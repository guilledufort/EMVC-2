#!/bin/bash

# number of parallel threads
parallel_threads=5

# dataset
dataset="$1"

# get the filename from the input file
filename=$dataset/files_list.txt

# function to download a link in parallel
download_link() {
  link="$1"

  # extract the filename from the link
  file=$(echo $link | awk -F/ '{print $NF}')

  # check if the file already exists in the current folder
  if [ ! -f "$dataset/$file" ] || [ ! -s "$dataset/$file" ]; then
    # if not, download the file using curl
    curl -o "$dataset/$file" $link
  else
    # if it does, print a message indicating that it already exists
    echo "$file already exists in the $dataset folder"
  fi
}

# loop through each line in the input file
while IFS= read -r line
do
  # download the link in a new background process
  download_link "$line" &

  # limit the number of parallel processes
  if (( $(jobs | wc -l) >= parallel_threads )); then
    wait -n
  fi
done < "$filename"

# wait for all remaining background processes to finish
wait

