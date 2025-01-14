# LR_on_Met
Codes for simulations and the real data application for the paper _Binary Regression and Classification with Covariates in Metric Spaces_.
* Below, '/home/linyn/' is the working directory.

## codes for simulation studies
* Running the following scripts IN SEQUENCE would produce all tables for simulation studies in the paper
* Before running, MAKE SURE the parameter 'main.path' in  every script ('simulation_manifold.R', 'compared_methods.R', 'stat_table.R' and './spdnet/task.py') is set to the working directory
* Note that, the first three commands are time-consuming. With our environment (20 CPU cores), it would take more than 12 DAYS to finish all the simulations

```
nohup Rscript /home/linyn/LR_on_Met/simulations/simulation_manifold.R null null 500 20 null	# calculating results for the proposed method
nohup Rscript /home/linyn/LR_on_Met/simulations/compared_methods.R null null 500 20 null	# calculating results for competing methods expect for SPDnet
nohup python /home/linyn/LR_on_Met/simulations/spdnet/task.py --case 0 --mu null --beta null	# calculating results for SPDnet
nohup Rscript /home/linyn/LR_on_Met/simulations/stat_table.R estimation 1 500			# generate Table 1 in the main paper and Table S1 in Supplementary Material
nohup Rscript /home/linyn/LR_on_Met/simulations/stat_table.R classification 1 500		# genearte Case 1 in Table 2 of the main paper and Tables S2-S3 of Supplementary Material
nohup Rscript /home/linyn/LR_on_Met/simulations/stat_table.R classification 2 500		# genearte Case 2 in Table 2 of the main paper and Case 2 in Tables S4 of Supplementary Material
nohup Rscript /home/linyn/LR_on_Met/simulations/stat_table.R classification 3 500		# genearte Case 3 in Table 2 of the main paper and Case 3 in Tables S4 of Supplementary Material
```

* './LR_on_Met/simulations/simulation_manifold.R': produce the results of the proposed method for the given settings, which would further be used to generate all tables in the paper. The arguments for the script include:
  - mu.types: string, types of the true Frechet mean (mu.star). It could be 'diag'(I_3), 'AR1'(\Sigma_{AR}) or 'null', where 'null' means iterating over both values.
  - beta.types: string, types of the true coefficient object (beta.star). It could be 'diag'(3I_3), 'AR1'(3\Sigma_{AR}) or 'null', where 'null' means iterating over both values.
  - R: positive integer, number of replicates.
  - ncore: positive integer, number of CPU cores to use.
  - case: positive integer, number of the model case. It could be 1, 2, 3, or null, where 'null' means iterating over all three values.
* For example, the following command would calculating the simulation results of the proposed method for Case 1 with mu.star=I_3 and beta.star=3\Sigma_{AR} 
```
nohup Rscript /home/linyn/LR_on_Met/simulations/simulation_manifold.R diag AR1 500 20 1	
```

* './LR_on_Met/simulations/compared_methods.R': produce the results of the competing methods expect for SPDnet for the given settings, which would further be used to generate all tables in the paper. The arguments for the script include: 
  - mu.types: string, types of the true Frechet mean (mu.star). It could be 'diag'(I_3), 'AR1'(\Sigma_{AR}) or 'null', where 'null' means iterating over both values. 
  - beta.types: string, types of the true coefficient object (beta.star). It could be 'diag'(3I_3), 'AR1'(3\Sigma_{AR}) or 'null', where 'null' means iterating over both values.
  - R: positive integer, number of replicates.
  - ncore: positive integer, number of CPU cores to use.
  - case: positive integer, number of the model case. It could be 1, 2, 3, or null, where 'null' means iterating over all three values.
* For example, the following command would calculating the simulation results of the competing methods expect for SPDnet for Case 1 with mu.star=I_3 and beta.star=3\Sigma_{AR} 
```
nohup Rscript /home/linyn/LR_on_Met/simulations/compared_methods.R diag AR1 500 20 1	
```

* './LR_on_Met/simulations/spdnet/task.py': produce the results of SPDnet for the given settings, which would further be used to generate all tables in the paper. The arguments for the script include: 
  - case: positive integer, number of the model case. It could be 1, 2, 3, or 0, where 0 means iterating over all three values.
  - mu: string, types of the true Frechet mean (mu.star). It could be 'diag'(I_3), 'AR1'(\Sigma_{AR}) or 'null', where 'null' means iterating over both values. 
  - beta: string, types of the true coefficient object (beta.star). It could be 'diag'(3I_3), 'AR1'(3\Sigma_{AR}) or 'null', where 'null' means iterating over both values.
* For example, the following command would calculating the simulation results of SPDnet for Case 1 with mu.star=I_3 and beta.star=3\Sigma_{AR} 
```
nohup python /home/linyn/LR_on_Met/simulations/spdnet/task.py --case 1 --mu diag --beta AR1
```

* './LR_on_Met/simulations/stat_table.R': Generate all tables in the paper. Note that this command depends on the results by 'simulation_manifold.R', 'compared_methods.R' and './spdnet/task.py'. The arguments for the script include: 
  - task: string, 'Estimation' or 'Classification' for estimation or classification performance.
  - case: positive integer, number of the model case. It could be 1, 2, 3.
  - R: positive integer, number of replicates.
* For example, the following command would generate Table 1 in the main paper and Table S1 in Supplementary Material
```
nohup Rscript /home/linyn/LR_on_Met/simulations/stat_table.R estimation 1 500
```


## codes for data application
* Running the following scripts IN SEQUENCE would produce all tables/figures for data application in the paper
* Before running, MAKE SURE the parameter 'main.path' in  every script ('real_analysis.R', 'compared_methods.R', 'nullity_test.R') is set to the working directory
* Note that, these commands are time-consuming. With our environment (20 CPU cores), it would take more than ONE WEEK to finish all the analysis
```
nohup Rscript /home/linyn/LR_on_Met/real_data/fMRI/real_analysis.R	# calculating results for the proposed method, and generating Figures 1-2 in the main paper
nohup Rscript /home/linyn/LR_on_Met/real_data/fMRI/nullity_test.R	# perform the permutation test for the nullity of the covariate effects, and generating Figure S1 in Supplementary Material
nohup python /home/linyn/LR_on_Met/real_data/fMRI/spdnet/task.py --data-dir /home/linyn/LR_on_Met/real_data/fMRI/matlab-data	# calculating results for SPDnet
nohup Rscript /home/linyn/LR_on_Met/real_data/fMRI/compared_methods.R	# calculating results for competing methods expect for SPDnet, and generating Table 3 in the main paper
```

