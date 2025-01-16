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
		const std::chrono::time_point<std::chrono::system_clock> start_time = std::chrono::system_clock::now();
	
	public:
		Benchmark(std::string method, std::filesystem::path reference_path, std::filesystem::path query_path);
		void write(int read_num);
};

#endif
