#Instructions
You need to extract the logs before running the notebook and also awk the files to get something that a computer can read

tar -xvf raw_53bits.tar.gz 

awk '/global/ {getline;{print}}' log.gpurun > series_all_notrunc
