#!/bin/bash

# Check if a directory is provided as an argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

# Define the path to your executable
executable="./raytracer"

# Define the directory containing input files
directory="$1"

# Check if the directory exists
if [ ! -d "$directory" ]; then
    echo "Error: Directory '$directory' does not exist."
    exit 1
fi

# Get a list of input files in the directory
inputs=("$directory"/*)

# Loop through the input files
for input in "${inputs[@]}"
do
    echo "Running $input..."
    
    # Measure the time taken to execute the command
    start_time=$(date +%s.%N)
    $executable "$input"
    end_time=$(date +%s.%N)
    
    # Calculate and print the elapsed time
    elapsed_time=$(echo "$end_time - $start_time" | bc)
    echo "Elapsed time: $elapsed_time seconds"
    echo "--------------------------------"
done

