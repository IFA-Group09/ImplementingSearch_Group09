#include <sstream>
#include <ranges>
#include <algorithm>
#include <tuple>
#include <unordered_set>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

#include "benchmark.hpp"

struct match_hash { 
  size_t operator()(const std::tuple<int, int, int> &val) const { 
    return std::get<0>(val) ^ std::get<1>(val) ^ std::get<2>(val); 
  } 
}; 

bool verify(std::vector<seqan3::dna5> const& ref, std::vector<seqan3::dna5> const& query, int start_position, int max_mismatches) {
    int mismatches = 0;
    for (auto j = 0; j < query.size(); j++) {
	if (mismatches > max_mismatches)
		break;

	if (ref[start_position+j] != query[j]) {
		mismatches++;
	}
    }

    return mismatches <= max_mismatches;
}

int main(int argc, char const* const* argv) {
    seqan3::argument_parser parser{"fmindex_pigeon_search", argc, argv, seqan3::update_notifications::off};

    parser.info.author = "SeqAn-Team";
    parser.info.version = "1.0.0";

    auto index_path = std::filesystem::path{};
    parser.add_option(index_path, '\0', "index", "path to the query file");

    auto reference_file = std::filesystem::path{};
    parser.add_option(reference_file, '\0', "reference", "path to the reference file");

    auto query_file = std::filesystem::path{};
    parser.add_option(query_file, '\0', "query", "path to the query file");

    auto number_of_queries = size_t{100};
    parser.add_option(number_of_queries, '\0', "query_ct", "number of query, if not enough queries, these will be duplicated");

    auto number_of_errors = uint8_t{0};
    parser.add_option(number_of_errors, '\0', "errors", "number of allowed hamming distance errors");

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

    // loading fm-index into memory
    using Index = decltype(seqan3::fm_index{std::vector<std::vector<seqan3::dna5>>{}}); // Some hack
    Index index; // construct fm-index
    {
        seqan3::debug_stream << "Loading 2FM-Index ... " << std::flush;
        std::ifstream is{index_path, std::ios::binary};
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
        seqan3::debug_stream << "done\n";
    }

    // duplicate input until its large enough
    while (queries.size() < number_of_queries) {
        auto old_count = queries.size();
        queries.resize(2 * old_count);
        std::copy_n(queries.begin(), old_count, queries.begin() + old_count);
    }
    queries.resize(number_of_queries); // will reduce the amount of searches

    seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}};

    auto benchmark = Benchmark("fmindex_pigeon", reference_file, query_file, 0);
    int read_num = 0;
    for (auto& query : queries) {
	std::unordered_set<std::tuple<int, int, int>, match_hash> match_results;
	int piece_size = query.size()/(number_of_errors+1);
	int first_offset = query.size() % (number_of_errors+1);
	std::vector<std::span<seqan3::dna5>> pieces;
	for (auto i = 0; i < (number_of_errors+1); i++) {
		int start;
		int end;
		if (i == 0) {
			start = i*piece_size;
			end=piece_size+first_offset;
		} else {
			start = (i*piece_size)+first_offset;
			end=piece_size;
		}
		auto piece = std::views::counted(query.begin()+start, end);
		pieces.push_back(piece);
	}
	auto results = seqan3::search(pieces, index, cfg);

    	for (auto && result : results) {
		match_results.insert(std::make_tuple(result.reference_begin_position(), result.query_id(), result.reference_id()));
	}

	for (auto& [match_position, piece_id, reference_id] : match_results) {
		// if we cannot have possibly found a match within reference bounds then skip
		if ( ((match_position - ((piece_id*piece_size)+first_offset)) < 0) || (match_position + (piece_id*piece_size)+first_offset > reference[reference_id].size())) {
			continue;
		}

		int matched_pieces = 0;
		for (auto i=0; i < piece_id; i++) {
			int offset = piece_size * (piece_id - i);
			if (i == 0)
				offset += first_offset;
			
			if (match_results.contains(std::make_tuple(match_position-offset, i, reference_id))) {
				matched_pieces++;
			} 
		}

		for (auto i=piece_id+1; i < pieces.size(); i++) {	
			if (match_results.contains(std::make_tuple(match_position+(piece_size*i), i, reference_id))) {
				matched_pieces++;
			}
		}

		if (matched_pieces == pieces.size()-1) {
			seqan3::debug_stream << "exact match for '" << query << "' at position: " << match_position << "\n";
		} else if (matched_pieces == pieces.size()-2) {
			//seqan3::debug_stream << "Verifying partial match\n";
			if (verify(reference[reference_id], query, (match_position-((piece_size*piece_id)+first_offset)), number_of_errors)) {
				seqan3::debug_stream << "partial match for '" << query << "' found at position: " << match_position-(piece_size*piece_id) << "\n";
			}
		}
	}

	if (read_num % 10 == 0) {
		benchmark.write(read_num);
	}
	read_num++;

    }


    return 0;
}
