# `CODEML` analysis

## 1. Pick rate prior

We will use a vague gamma distribution centered on a mean evolutionary rate estimated by considering the tree height (molecular distance in substitutions per site) and the mean age for the root of our phylogenies (time unit = 1000 Mya). As the rooted phylogeny inferred by [Moody et al. 2022](https://elifesciences.org/articles/66695) has information about the branch lengths, we can use [our R in-house script](scripts/calculate_rateprior.R) to calculate the corresponding tree heights. We also have a calibration to constrain the root age: a minimum age of 3,347 Ma and a maximum age of 4,520 Ma; which mean we will use as an approximate age for the root of our phylogeny to estimate the mean evolutionary rate.

By setting a vague shape ($\alpha=2$) for the gamma distribution that we will use as a rate prior, we can account for the uncertainty surrounding our mean rate estimate. If we had more knowledge on the mean rate, however, we should use a narrower prior with a larger $\alpha$ that better represents our prior information.

Now, we have all the information we need to calculate the $\beta$ parameter for the gamma distribution! We have written the [R script `calculate_rateprior.R`](scripts/calculate_rateprior.R) to carry out all the tasks mentioned above. You can open this file in RStudio to find out the value of $\beta$ and plot the final prior on the rates. A summary of what you will find in the script is described below:

```text
First, we know that the molecular distance (tree height, distance from the root to present time) is equal to the mean evolutionary rate (in substitutions per site per time unit) times the age of the divergence time at the root (in time unit, which we can define later). If we have estimated our phylogeny, and therefore have estimated the value of the branch lengths, we will be able to calculate the tree height. The units of the tree height will be the following:

tree_height = rate * root_age --> units_tree_height = subst/site/time * time = subst/site

There are various functions we can use to calculate the tree heigt. We have chosen the R function `phytools::nodeHeights`. The maximum height calculated by this function corresponds to the length from the root to the heighest tip.

After calculating the tree height of our phylogeny (in subst/site) and considering the age of the root based on the fossil evidence (time unit = 1 Ga = 100 Mya = 1e9 years), we can get a rough estimate for the mean rate:

Time unit = 1 Ga (mean root age in Ga) --> mean_rate = tree_height / root_age = (subst/site) / (Ga) = subst/site/Ga (time unit = 1 Ga = 1e9 years) 

We also know that the mean of the gamma distribution that we want to use as rate prior is our parameter of interest: the mean evolutionary rate. Therefore:

mean_Gamma = mean_rate = alpha / beta 
Time unit = 1 Ga: mean_rate = alpha / beta --> beta = alpha / mean_rate = 2 / mean_rate

The calibrated tree needs to incorporate the age constraints in the same time unit used to infer the mean evolutionary rate and establish the rate prior (i.e., do not forget to scale the calibrations accordingly if needed!). 
```

If you run the [R script `calculate_rateprior.R`](scripts/calculate_rateprior.R), you will see how all the steps described above take place and a new PDF file with the prior distribution to be used will be generated in a new directory called `out_RData`.

We have then updated out [template control file](control_files/prepcodeml.ctl) with the $\alpha$ and $\beta$ parameters for the gamma distribution as a rate prior (same distribution used regardless of the tree hypothesis being tested). Note that several options in this control file will be subsequently modified to fit the analysis with this dataset (i.e., you will see some options that have flags in capital letters, which will be replaced with the correct value for said option). Given how shallow this tree is, the clock may be expected to be seriously violated, and thus we have fixed a mean for the `sigma2` parameter (i.e., variation in the clock) as 0.1 using a gamma prior with $\alpha=1$ and $\beta=10$: `sigma2_gamma 1 10​`.

## 2. Set up the file structure

To approximate the likelihood calculation to speed up timetree inference ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) with `MCMCtree`, we shall need the branch lengths (vector), the gradient (vector), and the Hessian (matrix). We will use `CODEML` for this purpose as our dataset is an amino acid alignment. The file structure we will use is the following:

```text
anoxphoto/ # Our working directory
  |
  |- alignments/
  |    |- X/ # Directory for alignment X -- have as many directories as alignments or tree hypotheses
  |       
  |- control_files/ # Pre-defined control file with flags to be later replaced with specific settings
  |
  |- Hessian/
  |    |- X # Directory for alignment X -- have as many directories as alignments or tree hypotheses
  |          
  |- pipelines_Hessian # Directory where the pipeline to run `CODEML` will be executed
  |
  |- scripts # Scripts used to prepare control files to run `CODEML
  |
  |- trees
      |- calibrated   # Directory with the calibrated tree/s for `MCMCtree`
      |- uncalibrated # Directory with the uncalibrated tree/s for `CODEML`
```

To create the `anoxphoto` file structure, we can run the following commands on a local PC before transferring the full file structure to a HPC:

```sh
# Run the following commands from 
# directory `00_CODEML`
mkdir -p HPC/anoxphoto
cd HPC/anoxphoto
# One directory for each calibrated tree file
# to be used by `MCMCtree`
num_trees_calib=2
for i in `seq 1 $num_trees_calib`
do
mkdir -p trees/{uncalibrated,calibrated/$i}
done
# Only one uncalibrated tree (i.e., same tree
# topology, only calibrations change)
mkdir -p Hessian/1/prepare_codeml
mkdir -p alignments # Only one alignment
mkdir -p control_files
mkdir -p pipelines_Hessian
mkdir scripts
```

Once the file structure is created, we can now populate it with the input files we have generated some minutes ago: alignment file, tree files, and control file. We will also add the `lg.dat` file, which has the matrix to enable `CODEML` to use the LG substitution model. You can transfer the files to the HPC as you prefer, but below we show an example of how you could do this using `rsync`, a system that we will keep as an example throughout the rest of the tutorial:

```sh
# Run from `HPC/anoxphoto`
# Copy alignment
cp ../../../../00_data_formatting/01_inp_data/ortho57_aln.phy alignments/
# Now, transfer calibrated and uncalibrated trees
name_dat=( 'withLACA' 'withoutLACA' )
count=-1
num_trees=2
for i in `seq 1 $num_trees`
do
count=$(( count + 1 ))
cp ../../../../00_data_formatting/01_inp_data/anoxphoto_${name_dat[count]}_calib_MCMCtree.tree trees/calibrated/$i
done
# Transfer uncalibrated tree file
cp ../../../../00_data_formatting/01_inp_data/anoxphoto_uncalib.tree trees/uncalibrated
# Next, copy control files
cp ../../control_files/* control_files/
# Last, copy the in-house bash scripts with our pipeline
cp ../../scripts/*sh scripts/
# Once everything is ready, you can transfer this directory to your cluster!
# One way of doing this is by using `rsync`, but you may use other approaches
# Below, you will find an example of the `rsync` commands you should run once
# you replace the tags with your own credentials
# Move first one dir back so you are inside `HPC`, then run the `rsync` command
cd ../
rsync -avz --copy-links anoxphoto <uname>@<server>:<path_to_your_wd_in_HPC>
```

----

**IMPORTANT INFORMATION REGARDING THE PAML VERSION INSTALLED ON OUR HPC**
We have compiled the PAML programs available for the latest version of this software in our HPC (i.e., PAML v4.10.7) given that the cross-bracing approach is implemented in such a version. There are two ways that you could follow to use PAML software on your HPC:

* **If you have an older version exported to your path variable and do not want to change it**: firstly, make sure that you compile `MCMCtree` after modifying the source code: `NS` needs to be increased as there are more than 500 taxa in our alignment! To make these changes, please open the `mcmctree.c` file inside the `src` directory and change the line where `NS` is defined to read as it follows: `#define NS            2000`. Then, please save the changes and compile `MCMCtree` following the PAML installation guidelines given on [the PAML GitHub repository](https://github.com/abacus-gene/paml/wiki#installation). Once the software is compiled, please rename the `mcmctree` and `codeml` binaries to `mcmctree4.10.7` and `codeml4.10.7`, respectively, and save them in your `anoxphoto` working directory so that they are launched using relative paths.
* **If you want to export the latest PAML version to your path variable**: as explained above, you will need to compile `MCMCtree` after modifying the source code: `NS` needs to be increased as there are more than 500 taxa in our alignment! To make these changes, please open the `mcmctree.c` file inside the `src` directory and change the line where `NS` is defined to read as it follows: `#define NS            2000`. Then, please save the changes and compile `MCMCtree` following the PAML installation guidelines given on [the PAML GitHub repository](https://github.com/abacus-gene/paml/wiki#installation). Then, please export the path to the `bin` directory where you should have all the PAML binaries after compilation took place. You will be able to launch these programs by typing `mcmctree` and `codeml` from any location in your file structure if the path has been properly exported.

Please note that `MCMCtree` and `CODEML` are the PAML programs that we will use for timetree inference. In this step-by-step tutorial, you shall see how we call these programs via relative paths (if you chose the former) or absolute paths (if you chose the latter). Given that `conda` does not yet have the latest PAML version available, these are the two possibilities you have to work with the latest PAML version on your HPC.

----

Now, we need to generate other input files to estimate the branch lengths, the Hessian, and the gradient: the input control files for `CODEML`. To do this in a reproducible manner, you can use the [script `generate_prepcodeml.sh`](scripts/generate_prepcodeml.sh), which you can find in the [`01_PAML/00_CODEML/scripts`](01_PAML/00_CODEML/scripts) and which you should have just transferred to your HPC. Now, connect to your server and run the next code snippet, where you will execute this script. Specifically, the [`generate_prepcodeml.sh` script](scripts/generate_prepcodeml.sh) needs one argument: the amount of alignment files to be analysed. As we only have one alignment file, we will use `1` as the argument:

```sh
# Run from `anoxphoto/scripts` in the HPC
# Please change directories until
# you are there
# Then, run the following commands
chmod 775 *sh
# In this case, there is only 1 tree topology
num_tt=1
for i in `seq 1 $num_tt`
do
./generate_prepcodeml.sh $i
done
```

To make sure that all the paths have been properly extracted, you can run the following code snippet:

```sh
# Run from `anoxphoto/Hessian` dir on your local PC
# Please change directories until you are there
# Then, run the following commands
grep 'seqfile' */prepare_codeml/*ctl
grep 'treefile' */prepare_codeml/*ctl
grep 'aaRatefile' */prepare_codeml/*ctl
```

## 3. Run `CODEML`

### Preparing input files

Now that we have the input files (alignment and tree files) and the instructions to run `CODEML` (control file) in our HPC, we will be manually running `MCMCtree` inside each `prepare_codeml` directory (see file structure above) in a special mode that launches `CODEML` for the sole purpose we want: to infer the vectors and matrix required to approximate the likelihood calculation.

```sh
# Run `MCMCtree` from
# `anoxphoto/Hessian/1/prepare_codeml`
# dir on the HPC
# Please change directories until
# you are in there
# The first command to change directories 
# will work if you are still in 
# `main/Hessian`, otherwise ignore and 
# move to such directory with the command
# that best suits your current directory
cd 1/prepare_codeml
##>> RUN THE NEXT COMMAND IF USING LOCAL BINARY
##>> NOT EXPORTED TO THE PATH'S SYSTEM
../../../mcmctree4.10.7 prepcodeml.ctl
##>> END
##>> RUN THE NEXT COMMAND IF PROGRAM EXPORTED
##>> TO THE PATH'S SYSTEM
mcmctree prepcodeml.ctl
##>> END
```

Firstly, you will see that `MCMCtree` starts parsing the first locus. Then, you will see something like the following printed on your screen (some sections may change depending on the PAML version you have installed on your cluster!):

```text
*** Locus 1 ***
running codeml tmp0001.ctl

AAML in paml version 4.9h, March 2018
ns = 700        ls = 8144
Reading sequences, sequential format..
Reading seq #700: GCA_900620265       
Sequences read..

8144 site patterns read, 8152 sites
Counting frequencies..

  1957200 bytes for distance
  2606080 bytes for conP
   260608 bytes for fhK
  5000000 bytes for space
```

> [!NOTE]
> Even though you see `version 4.9h` written in the example above, we are using this strategy to generate the input files, not to run `CODEML` itself! Therefore, it does not matter if you run either older versions or the latest version of this program!

As soon as you see the last line, you will see that various `tmp000X*` files will have been created, and hence you can stop this run by typing `ctrl+C` on the terminal that you have used to run the command. Once you have done this, you can check that the control file you will later need has been created:

```sh
# Run from the `anoxphoto/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */*/tmp0001.ctl | wc -l # You should get as many hypotheses as you have, in this case 1
```

Please note that, when we ran the commands above, we were not interested in running `CODEML` or `MCMCtree`. We just wanted to execute `MCMCtree` with option `usedata = 3` so that it generates the `tmp000*` files that `CODEML` will later need to estimate the branch lengths, the gradient, and the Hessian. We do this analysis in two steps given that there are restrictions in the HPC we are using that do not allow us to run `CODEML` + `MCMCtree` in one unique job within a reasonable amount of time. In addition, we want to modify some settings in the control file that is automatically generated when enabling `usedata = 3` so that they match what we want to do for our inference. In a nutshell, this is what you will be doing:

1. Run `MCMCtree` to generate the `tmp000*` files.
2. Modify the `tmp0001.ctl` file according to the settings we want to enable to analyse our dataset with `CODEML`.
3. Run `CODEML` using the `tmp000*` files so that it branch lengths, the gradient, and the Hessian can be estimated and saved in a file called `rst2`.
4. Generate the final `in.BV` file for our dataset, which will be later used by `MCMCtree` to approximate the likelihood calculation.

Once all `tmp000*` files are generated for all alignments, we need to make sure that the following options have been enabled:

1. The (absolute or relative) path to the `lg.dat` file with the AA substitution model should be the argument of option `aaRatefile` in the `tmp000*.ctl`.
2. Following the `Tutorial 4: Approximate likleihood with protein data`, one of the sections in the [`MCMCtree Tutorials` document](http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf), we set the template control file to use gamma rates among sites instead of the default model, which uses no gamma rates. As we had already defined these options in the template control file, these are already enabled in the `tmp000*.ctl` file generated after running `MCMCtree`. In other words, the control file should already have (i) `fix_alpha = 0` and `alpha = 0.5` to enable the estimation of alpha for the gamma distribution for variable substitution rates across sites by starting the search of the value of $\alpha$ at 0.5 and (ii) `ncatG = 4` with the number of categories for this distribution equal to 4. Nevertheless, we can always double check that this is indeed the case to make sure that everything is fine.
3. According to the settings above, we will be using the "LG+F+G4" empirical model for amino acid data, which is enabled by setting `model = 3`, `ncatG = 4`, `aaRatefile = <path_to_lg_matrix>/lg.dat`.
4. In addition, you need to make sure that option `method = 1` is enabled, which will speed up the calculations.

We can run the next code snippet to very that the four requirements aforementioned are met:

```sh
# Run from the `anoxphoto/Hessian` dir on your local PC
# Please change directories until you are there
# Then, run the following commands
sed -i 's/method\ \=\ 0/method\ \=\ 1/' */*/tmp0001.ctl
grep 'method = 1' */*/tmp0001.ctl | wc -l # You should get as many as tree hypotheses as you have
grep 'alpha' */*/tmp0001.ctl   # You should see `fix_alpha = 0` and `alpha = 0.5`
grep 'lg.dat' */*/tmp0001.ctl  # You should see the absolute path to the `lg.dat` file in your system
grep 'ncatG' */*/tmp0001.ctl   # You should see `ncatG = 4`
grep 'model' */*/tmp0001.ctl   # You should see `model = 3` (i.e., empirical+F model)
```

### Executing `CODEML`

We can now run `CODEML` given that we have the control file ready as well as all the required input files!

We have created a template bash script with flags (i.e., see script `pipeline_Hessian_CODEML_template.sh` in the [`scripts` directory](01_PAML/00_Hessian/scripts/pipeline_Hessian_CODEML_template.sh)), which will be replaced with the appropriate values by another bash script (`generate_job_CODEML.sh`, also saved in the [`scripts` directory](01_PAML/00_Hessian/scripts/generate_job_CODEML.sh)). Please note that the second bash script will edit the template bash script according to the data alignment/s that will be analysed. We had already transferred these scripts to the HPC when setting up our file structure. Therefore, we just need to execute the following code snippet there:

```sh
# Run from `anoxphoto` dir on your HPC
# Please change directories until you are there
# Then, run the following commands
home_dir=$( pwd )
cd scripts
chmod 775 *sh
num_trees=1
##>> RUN THE NEXT COMMAND ONLY IF SGE-BASED SYSTEM
# Arg1: Number of alignments
# Arg2: Path to the pipeline directory
# Arg3: Name of the working directory (i.e., `anoxphoto` in this analysis)
# Arg4: Name of the executable file for CODEML. E.g., `codeml4.10.7`, `codeml`, etc.
# Arg5: Boolean, PAML exported to the path? `Y` of `N`
# Arg6: Requested RAM. E.g., `2G`
##>> END IF SGE-BASED SYSTEM
./generate_job_CODEML.sh $num_trees $home_dir/pipelines_Hessian anoxphoto codeml4.10.7 N "20G"
##>> RUN THE NEXT COMMAND ONLY IF SLURM-BASED SYSTEM
# Arg1: Number of alignments
# Arg2: Path to the pipeline directory
# Arg3: Name of the working directory (i.e., `anoxphoto` in this analysis)
# Arg4: Name of the account to use, if any. If not, just enter NA
# Arg5: Name of the partition to use, if any. If no partition, enter NA
# Arg6: Name that will be used to tag this job
# Arg7: Name of the executable file for CODEML. E.g., `codeml4.10.7`, `codeml`, etc.
# Arg8: Boolean, PAML exported to the path? `Y` of `N`
# Arg9: Requested RAM. E.g., `2G`
# Arg10: Ilimited time? `Y` or `N`
./generate_job_CODEML_SLURM.sh $num_trees $home_dir/pipelines_Hessian anoxphoto GELY012201 taw abce_codeml codeml Y "20G" Y
##>> END IF SGE-BASED SYSTEM
```

Next, we will go to the `pipelines_Hessian` directory and run the script that will have been generated using the commands above:

```sh
# Run from `anoxphoto/pipelines_Hessian` dir on your HPC
# Please change directories until you are there
# Then, run the following commands
#
# If you list the content of this directory,
# you will see the pipeline you will need 
# to execute in a bash script called
# `pipeline_Hessian.sh`
ll *
# Now, execute this bash script
chmod 775 *sh
qsub pipeline_Hessian.sh
```

Once `CODEML` finishes, we are ready to generate the `in.BV` file that we will later use when running `MCMCtree` to approximate the likelihood calculation:

```sh
# Run from dir `anoxphoto/Hessian/` dir on your HPC
# Please change directories until you are there
# Then, run the following commands
printf "\nGenerating in.BV files for dir 1 ... ...\n\n"
cp 1/rst2 in.BV
```

Next, we can transfer the output generated by `CODEML` to our HPC so that we can keep a backup:

```sh
# Run from `00_CODEML_conc` in your PC
mkdir out_CODEML
cd out_CODEML
rsync -avz --copy-links <uname>@<server>:<path_to_your_wd_in_HPC>/anoxphoto/Hessian .
rsync -avz --copy-links <uname>@<server>:<path_to_your_wd_in_HPC>/anoxphoto/pipelines_Hessian .
# Remove unnecessary empty output files
rm pipelines_Hessian/*sh.o*
```

We can now proceed to timetree inference with `MCMCtree`. [You can click this link to move to the next `README.md` file](../01_MCMCtree/README.md)!
