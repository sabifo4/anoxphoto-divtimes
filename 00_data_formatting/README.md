# Generate input data

Before proceeding with timetree inference, we need to make sure that:

1. The alignment files are in PHYLIP format and easy to read (i.e., ideally one sequence per line).
2. The tree file is in Newick format.

## Alignment files

The first step is to generate the input sequence alignment! In the [`00_raw_data`](00_raw_data/alignment/ind_genes), you will find a directory called [`ind_genes`](00_raw_data/alignment/ind_genes), where the individual gene alignments that were generated by [Moody et al. 2022](https://elifesciences.org/articles/66695) can be found. We will run the following code snippet to generate a concatenated alignment in PHYLIP format:

```sh
# Run from `00_raw_data/alignment`
mkdir ind_to_conc
cd ind_to_conc
maln_dir=$( pwd )
cd ../ind_genes
home_dir=$( pwd )
# 1. Move to `src`, where we shall uncompress the file
#    with the `fasta-phylip-partitions` to generate the
#    concatenated alignment file.
#    For more details about how to use this pipeline, you
#    can find this in the following link:
#    https://github.com/sabifo4/fasta-phylip-partitions/blob/main/README.md
cd ../../../../src
tar -xvf fasta-phylip-partitions.tar.gz
chmod 775 fasta-phylip-partitions/src/*sh
chmod 775 fasta-phylip-partitions/src/Tools/*
cd fasta-phylip-partitions/src/
base_dir=$( pwd )
cd $home_dir
# 2. Copy one-line FASTA files in a unique directory individually
#    to generate individual PHYLIP alignment as well as in a
#    separate directory to then generate the main alignment
#    following PHYL format
for i in *BMGE.faa
do
name=$( echo $i | sed 's/\_mafft..*//' )
cp $i $maln_dir/$name.fasta
done
# 3. Create `species_names.txt`. In order to create a file listing the species 
#    in the alignment, we will use the `57_BMGE_mafft_concat_SR4.fasta`
cd $home_dir
grep '>' ../nuc_conc/57_BMGE_mafft_concat_SR4.fasta | sed 's/>//g' > $maln_dir/species_names.txt
# Now, run the `fasta-phylip-partitions` pipeline!
# In essence, the first argument is the current 
# directory ("."), the second argument the tag ID for the job
# ("cytb"), and then a tag specifying if the alignment 
# needs to be partitioned into CPs, which we do not want,
# and hence use "partN".
# If you do not have this pipeline on your system, you can
# run the code below as it is using the source code that has already
# been added in the `src` directory.
# NOTE: If you are running this code in a directory that you have
# synched to Dropbox or another cloud storage and you have issues,
# just move the folder out of the synched directory and run the 
# command below again
cd ../ind_to_conc
$base_dir/Run_tasks.sh . ortho57 partN
# 4. Generate a directory for input data
cd $maln_dir/../../../
mkdir 01_inp_data
cp 00_raw_data/alignment/ind_to_conc/phylip_format/02_concatenated_alignments/ortho57_concat.aln 01_inp_data/ortho57_aln.phy
```

Now that we have formatted the alignment file, we are ready to prepare the tree files!

## Tree files

We are going to use the main 700-taxa tree with all the taxonomical information and integrate both geological and fossil record evidence to constrain some of the node ages of the 700-taxa phylogeny. We will be using two tree files generated by [Moody et al. 2022](https://elifesciences.org/articles/66695) with the constrained species tree so that we can generate the PHYLIP tree files in the correct format for the subsequent PAML analyses:

```sh
# Run from `00_raw_data/trees`, a directory that will be created
# by the R script aforementioned
cp constrained_sp_tt/constrained_species_rooted.tre tree_bl.tree
name=$( echo *bl.tree | sed 's/\_bl..*//' )
sed 's/\:[0-9]\.[0-9]*//g' *bl.tree | sed 's/\:[0-9]//g' | sed 's/\:[0-9]*\.[0-9]*e-[0-9]*//g' | sed 's/\:[0-9]e-[0-9]*//g' | sed 's/E-[0-9]*//g' > $name"_nobl.tree"
sed -i '1 i\700 1' $name"_nobl.tree"
# Now, copy the tree with all taxonomical info to manually include 
# annotations in `FigTree`
cp constrained_sp_tt/constrained_species_rooted.tre.renamed tree_detailed_taxonomy_annot.tree
# We have manually opened the file above in `FigTree` and included
# annotations to visualise which nodes in the tree are to be 
# calibrated. The resulting annotated tree has been saved in 
# various formats in the following tree files:
# --> `anoxphoto_annotate_FigTree`
# --> `anoxphoto_annotate_FigTree.nexus`
# --> `anoxphoto_annotat_FigTree.svg`
```

Now, we can use our [annotated phylogeny](00_raw_data/trees/anoxphoto_annotate_FigTree.svg) as a guideline to build [our calibration file](00_raw_data/calibs/Calibs_anoxphoto_withArchExclDPANN.txt). This file is formatted as follows:

* Header with eight columns:
  * name: tag that users want to use to refer to a specific node being calibrated. Spaces should not be used!
  * leaf1: one of the taxa that is used to identify the clade which MRCA corresponds to the node being calibrated. Users must use the tip label they have used to identify this taxon in their tree topology, no spaces!
  * leaf2: the second taxa that is part of the clade which MRCA corresponds to the node being calibrated.
  * MCMCtree: `MCMCtree` notation used to calibrate the MRCA node specified with `leaf1` and `leaf2` options aforementioned. Please write the calibration without spaces and within single quotation marks. E.g., for a soft-bound calibration with hard bounds, you could use `'B(3.225,4.520,1e-300,1e-300)'`.

> **EXAMPLE 1**: row in our text file to calibrate node "LUCA"

  ```text
  LUCA;GCA_000008085;GCA_000021645;'B(3.347,4.520,1e-300,1e-300)'
  ```

We will generate a copy of this file in which we will remove the row for the Crown-Archaea excluding DPANN calibration (i.e., `CG-ARCH-EXCLDPANN` label) to test the impact of including such calibration during timetree inference (see the resulting file [`Calibs_anoxphoto.txt`](00_raw_data/calibs/Calibs_anoxphoto.txt)).

Once both text files are ready, we can use [our in-house R script](scripts/Include_calibrations.R) to incorporate the `MCMCtree` formatted calibrations in the input tree files we will use during timetree inference analysis. Please note that time unit is 1Ga = 1000Ma. The resulting calibrated tree files will be automatically saved inside `01_inp_data` when you run [the `Include_calibrations` R script](scripts/Include_calibrations.R).

Now, we have the input alignment and tree files, so we can [get started with PAML software: let's run `CODEML` for branch lengths, gradient, and Hessian inference!](../01_PAML/00_CODEML/README.md)!
