#include "benchmark.hpp"

#include <fmindex-collection/fmindex-collection.h>
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

std::tuple<int, int> naive_binary_search(std::vector<seqan3::dna5>* query, std::vector<seqan3::dna5>* reference, std::vector<long unsigned int>* sa) {
	unsigned long int min_index = 0;
	unsigned long int max_index = sa->size();

	while (min_index < max_index) {
		auto c = (min_index + max_index)/2;
		auto ref_view = *reference | std::views::drop(sa->at(c)) | std::views::take(query->size()) | seqan3::views::to_char;
		auto query_view = *query | seqan3::views::to_char;
		if (std::ranges::lexicographical_compare(ref_view, query_view)) {
			min_index = c + 1;
		} else {
			max_index = c;
		}
	}

	auto first = min_index;
	max_index = sa->size();

	while (min_index < max_index) {
		auto c = (min_index + max_index)/2;
		auto ref_view = *reference | std::views::drop(sa->at(c)) | std::views::take(query->size()) | seqan3::views::to_char;
		auto query_view = *query | seqan3::views::to_char;

		if (std::ranges::lexicographical_compare(query_view, ref_view)) {
			max_index = c;
		} else {
			min_index = c + 1;
		}
	}
	auto last = max_index-1;
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

    auto quiet = false;
    parser.add_option(quiet, '\0', "quiet", "do not print matches");

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

    auto suffixarray = fmindex_collection::createSA64(std::span{reinterpret_cast<uint8_t const*>(reference.data()), reference.size()}, 1);
    int read_num = 0;
    auto benchmark = Benchmark("sa", reference_file, query_file, 0);
    for (auto& q : queries) {
        //!TODO !ImplementMe apply binary search and find q  in reference using binary search on `suffixarray`
        // You can choose if you want to use binary search based on "naive approach", "mlr-trick", "lcp"
	auto results = naive_binary_search(&q, &reference, &suffixarray);
	if (std::get<0>(results) >= 0) {
		for (auto i = 0; i < std::get<1>(results)-std::get<0>(results)+1; i++) {
			if (!quiet)
				seqan3::debug_stream  << q << "\n";
		}
	}

	if (read_num % 10 == 0) {
		benchmark.write(read_num);
	}
	read_num++;
    }

    return 0;
}
