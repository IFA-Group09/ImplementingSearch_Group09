#include <divsufsort.h>
#include <iostream>
#include <tuple>
#include <sstream>
#include <span>
#include <seqan3/utility/views/slice.hpp>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/alphabet/views/char_to.hpp>

std::tuple<int, int> naive_binary_search(std::vector<seqan3::dna5>* query, std::vector<seqan3::dna5>* reference, std::vector<saidx_t>* sa) {
	unsigned long int min_index = 0;
	unsigned long int max_index = sa->size();

	while (min_index < max_index) {
		auto c = (min_index + max_index) / 2;

		auto ref_view = *reference | std::views::drop(sa->at(c)) | seqan3::views::to_char;
		auto query_view = *query | seqan3::views::to_char;
		if (std::ranges::lexicographical_compare(ref_view, query_view)) {
			min_index = c + 1;
		} else {
			max_index = c;
		}
		
	}

	auto first = min_index;
	min_index = 0;
	max_index = sa->size();
	while (min_index < max_index) {
		auto c = (min_index + max_index) / 2;
		auto ref_view = *reference | std::views::drop(sa->at(c)) | seqan3::views::to_char;
		auto query_view = *query | seqan3::views::to_char;
		if (std::ranges::lexicographical_compare(query_view, ref_view)) {
			max_index = c;
		} else {
			min_index = c + 1;
		}
	}
	auto last = max_index;
	if ((first > last) || !(std::equal(reference->begin()+sa->at(first), reference->begin()+sa->at(first)+query->size(), query->begin()))) {
		return std::make_tuple(-1, -1);
	}

	return std::make_tuple(first, last);
}

int main(int argc, char const* const* argv) {
    seqan3::argument_parser parser{"suffixarray_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    auto number_of_queries = size_t{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "number of query, if not enough queries, these will be duplicated");

    try {
         parser.parse();
    } catch (seqan3::argument_parser_error const& ext) {
        seqan3::debug_stream << "Parsing error. " << ext.what() << "\n";
        return EXIT_FAILURE;
    }

    // loading our files
    auto reference_stream = seqan3::sequence_file_input{reference_file};
    auto query_stream     = seqan3::sequence_file_input{query_file};

    // read reference into memory
    // Attention: we are concatenating all sequences into one big combined sequence
    //            this is done to simplify the implementation of suffix_arrays
    std::vector<seqan3::dna5> reference;
    for (auto& record : reference_stream) {
        auto r = record.sequence();
        reference.insert(reference.end(), r.begin(), r.end());
    }

    // read query into memory
    std::vector<std::vector<seqan3::dna5>> queries;
    for (auto& record : query_stream) {
        queries.push_back(record.sequence());
    }

    // duplicate input until its large enough
    while (queries.size() < number_of_queries) {
        auto old_count = queries.size();
        queries.resize(2 * old_count);
        std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
    }
    queries.resize(number_of_queries); // will reduce the amount of searches

    // Array that should hold the future suffix array
    std::vector<saidx_t> suffixarray;
    suffixarray.resize(reference.size()); // resizing the array, so it can hold the complete SA

    //!TODO !ImplementMe implement suffix array sort
    //Hint, if can use libdivsufsort (already integrated in this repo)
    //      https://github.com/y-256/libdivsufsort
    //      To make the `reference` compatible with libdivsufsort you can simply
    //      cast it by calling:
    //      `sauchar_t const* str = reinterpret_cast<sauchar_t const*>(reference.data());`
    sauchar_t const* str = reinterpret_cast<sauchar_t const*>(reference.data());

    divsufsort((unsigned char *)str, &suffixarray[0], reference.size());

    // setup benchmarking output file
    std::ifstream benchmark_in;
    std::ofstream benchmark_f;
    benchmark_f.open("cpp_benchmark.csv", std::ios_base::app);

    benchmark_in.open("cpp_benchmark.csv");
    if (benchmark_in.peek() == std::ifstream::traits_type::eof()) {
	benchmark_f << "method,reads_file,time,read_n\n";
    }
    const auto start_time = std::chrono::system_clock::now();
    int read_num = 0;
    for (auto& q : queries) {
        //!TODO !ImplementMe apply binary search and find q  in reference using binary search on `suffixarray`
        // You can choose if you want to use binary search based on "naive approach", "mlr-trick", "lcp"
	auto results = naive_binary_search(&q, &reference, &suffixarray);
	if (std::get<0>(results) >= 0) {
		for (auto i = 0; i < std::get<1>(results)-std::get<0>(results)+1; i++) {
			seqan3::debug_stream  << q << "\n";
		}
	}

	if (read_num % 10 == 0) {
		benchmark_f << "sa," << query_file << "," << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now() - start_time).count() << "," << read_num << std::endl;
	}
	read_num++;
    }

    return 0;
}
