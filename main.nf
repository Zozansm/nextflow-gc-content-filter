#!/usr/bin/env nextflow

// Define input parameters
params.inputFile = null   // Path to the input FASTA file
params.cutoff = 0.5       // Default GC content threshold (fraction)

// Define the process to filter sequences by GC content
process filterGCContent {
    input:
    path fastaFile  // Input FASTA file

    output:
    path 'output.txt'  // Emit the result as output.txt

    script:
    """
    python3 -c "
from Bio import SeqIO

input_file = '${fastaFile}'
cutoff = ${params.cutoff}

with open(input_file, 'r') as handle, open('output.txt', 'w') as output:
    for record in SeqIO.parse(handle, 'fasta'):
        gc_content = (record.seq.count('G') + record.seq.count('C')) / len(record.seq)
        if gc_content > cutoff:
            output.write(f'>{record.id}\\n{record.seq}\\n')
    "
    """
}

// Define the workflow logic
workflow {
    // Ensure the input parameters are provided
    if (!params.inputFile) {
        error "Parameter '--inputFile' is required."
    }
    if (params.cutoff == null) {
        error "Parameter '--cutoff' is required."
    }

    // Create a path channel for the input FASTA file
    inputChannel = Channel.fromPath(params.inputFile)

    // Call the filtering process
    outputChannel = filterGCContent(inputChannel)

    // Display output file location
    outputChannel.view { "Filtered sequences written to: ${it}" }
}

