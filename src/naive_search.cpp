#include "benchmark.hpp"

#include <sstream>
#include <fstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

// prints out all occurences of query inside of ref
void findOccurences(std::vector<seqan3::dna5> const& ref, std::vector<seqan3::dna5> const& query, bool quiet) {
    for (long unsigned int i = 0; i <= ref.size()-query.size(); i++) {
	    for (long unsigned int j = 0; j <= query.size(); j++) {
		if (ref[i+j] != query[j])
			break;

		if ((j == query.size()-1) && !(quiet))
			seqan3::debug_stream << query << "\n";
	    }
    }
}

int main(int argc, char const* const* argv) {
    seqan3::argument_parser parser{"naive_search", argc, argv, seqan3::update_notifications::off};

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
    std::vector<std::vector<seqan3::dna5>> reference;
    for (auto& record : reference_stream) {
        reference.push_back(record.sequence());
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

    auto benchmark = Benchmark("naive", reference_file, query_file, 0);
    //! search for all occurences of queries inside of reference
    for (auto& r : reference) {
	int read_num = 0;
        for (auto& q : queries) {
            findOccurences(r, q, quiet);
	    if (read_num % 10 == 0) {
		    benchmark.write(read_num);
	    }
	    read_num++;

        }
    }

    return 0;
}
