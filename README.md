# Motif Mark

## The Program

This program will take a fasta file or a list of fasta files, as well as a list of motifs, and mark the motifs on each sequence in an svg plot for easy visualization of where the motifs fall on the sequence.  

### Requirements/Dependencies

This program requires python version 3.8.6, matplotlib version 3.3.4, and pycairo version 1.20.0.  It is recomended that these dependencies be installed in a conda environment.  This can be done through the following commands

`conda create env --name pycairo`
`conda activate pycairo`
`conda install -c conda-forge matplotlib`
`conda install -c conda-forge pycairo`

Or, the recomended alternative creating a conda environment from the "environment.yml" file included in this github.  This will create a conda environment called pycairo containing all required dependencies. 

```conda create --file environment.yml```

Following the creation of the conda environment, the environment must be activated prior to running the script.

```conda activate pycairo```

The program can then be run through the command line using the following command ./motif-mark.py -f fasta_file -m motif_file, -c "color_palette"

### Inputs

This program can take up to three input flags following the ./motif-mark.py command. -f (--files) and -m (--motifs) are both required.  --files must be followed by either a fasta file or a list of fasta files. The script will generate one image per fasta file and be named based on the prefix of that file name. The sequences contained within the fasta files must have proper fasta format. In addition the exon must be capitalized while the intron before and after the exon must be lowercase.  Currently, motif-mark only supports sequences in the intron, exon, intron format. The --motifs flag requires a singular text file containing a list of motifs (one per line) that have IUPAC naming conventions. Ambiguous nucletodies are allowed.  The third flag, -c (--color) is an optional parameter that allows you to change the matplotlib color palette being used.  Note the string put here must be a matplotlib color palette, otherwise this program will not function.  The default color palette is tab10.

## Output

The program will return one image per fasta file entered.  This image will contain a legend containing the colors associated with each motif.  Each sequence in the fasta file will have its own subplot labeled with the header line of the fasta file.  Each sequence will be proportional in length to the others and to the motifs.

If there are questions or concerns please contact the owner of this repository.