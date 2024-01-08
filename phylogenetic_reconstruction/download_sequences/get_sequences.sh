# Create directories for loci: 16s, 18s, 28s, COI, H3
mkdir 16s 18s 28s COI H3

# Loop through each folder in the current directory
for folder in $(ls -d */); do
    cd "$folder"
    IFS=$'\n'

    # Remove any existing download_sequences.fas and create an empty file
    rm -rf download_sequences.fas
    touch download_sequences.fas

    # Loop through each code in the codes.txt file within the folder
    for code in $(cat codes.txt); do
        # Download sequences using efetch and append them to download_sequences.fas
        esearch -db nucleotide -query "$code" | efetch -format fasta | grep ">" | sed -E "s/^>\w+.\w\s+(\w+\s+\w+)+.+/\>\1/g" > head.txt
        esearch -db nucleotide -query "$code" | efetch -format fasta | grep -v ">" > body.txt
        cat head.txt body.txt >> download_sequences.fas
        rm -rf head.txt body.txt
    done

    cd ..
done

# After obtaining all sequences for each locus, search for polymorphic species with genetic information to add to the phylogeny

# Loop through each folder in the current directory
for folder in $(ls -d */); do
    IFS=$'\n'
    cd "$folder"

    # Remove any existing head.txt
    rm -rf head.txt

    # Loop through each species in polymorphic_species.txt file in the parent directory
    for spe in $(cat ../polymorphic_species.txt); do
        rm -rf head.txt

        # Check if the species is not already in the current locus
        exists=$(grep "$spe" "$folder".fas | wc -l)
        if [[ $exists -eq 0 ]]; then
            # If not, download the sequence and add it to the locus
            esearch -db nucleotide -query "${spe}[Organism] ${folder}" | efetch -format fasta | grep ">" | sed -E "s/^>\w+.\w\s+(\w+\s+\w+)+.+/\>\1/g" > head.txt
            if [[ $(cat head.txt | wc -l) -gt 0 ]]; then
                esearch -db nucleotide -query "${spe}[Organism] ${folder}" | efetch -format fasta > "$spe"_sequences.fas
            else
                continue
            fi
        else
            continue
        fi
    done

    cd ..
done
