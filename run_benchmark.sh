
# First collect memory usage

READS_FILE=./data/illumina_reads_100.fasta.gz
REFERENCE_FILE=./data/hg38_partial.fasta.gz
FMINDEX_FILE=test.index

MEM_OUTPUT_FILE=cpp_memory_benchmark.csv
echo "method,reference_file,reads_file,read_num,mem_usage" > $MEM_OUTPUT_FILE
## Naive search
read_nums=(1000)
for read_num in "${read_nums[@]}"
do
	mem_usage=`\time -f "%M"  ./build/bin/naive_search --query $READS_FILE --reference $REFERENCE_FILE --query_ct $read_num 2>&1 | tail -n 1`
	echo "naive,${REFERENCE_FILE},${READS_FILE},${read_num},${mem_usage}" >> $MEM_OUTPUT_FILE
done



## SA search
read_nums=(1000 10000 100000 1000000)
for read_num in "${read_nums[@]}"
do
	mem_usage=`\time -f "%M" ./build/bin/suffixarray_search --query $READS_FILE --reference $REFERENCE_FILE --query_ct $read_num  2>&1 | tail -n 1`
	echo "sa,${REFERENCE_FILE},${READS_FILE},${read_num},${mem_usage}" >> $MEM_OUTPUT_FILE
done

## FM search
read_nums=(1000 10000 100000 1000000)
for read_num in "${read_nums[@]}"
do
	mem_usage=`\time -f "%M" ./build/bin/fmindex_search --query $READS_FILE --index $FMINDEX_FILE --query_ct $read_num 2>&1 | tail -n 1`
	echo "fm,${FMINDEX_FILE},${READS_FILE},${read_num},${mem_usage}" >> $MEM_OUTPUT_FILE
done

# remove any existing benchmark file
rm cpp_benchmark.cpp

# Now collect processing per read num
./build/bin/naive_search --query $READS_FILE --reference $REFERENCE_FILE --query_ct 1011 
./build/bin/suffixarray_search --query $READS_FILE --reference $REFERENCE_FILE --query_ct 1000011
./build/bin/fmindex_search --query $READS_FILE --index $FMINDEX_FILE --query_ct 1000011
