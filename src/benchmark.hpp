#ifndef BENCHMARK_HPP
#define BENCHMARK_HPP

#include <string>
#include <fstream>
#include <filesystem>

class Benchmark {
	private:
		std::string method;
		std::filesystem::path reference_path;
		std::filesystem::path query_path;
		std::ofstream benchmark_out;
		const std::chrono::high_resolution_clock::time_point start_time = std::chrono::high_resolution_clock::now();
		int number_of_errors;
	
	public:
		Benchmark(std::string method, std::filesystem::path reference_path, std::filesystem::path query_path, int number_of_errors);
		void write(int read_num);
};

#endif
