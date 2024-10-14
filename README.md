# Fast multi-trait analysis
The aim of the project is to introduce dimension reduction of outcome variables in multi-trait analysis.

# Data
## Scots pine field trial
The Scots pine example data and pedigree are provided in the directory data. The data contains all pre-adjusted phenotypes following Dutkowski and colleagues [[1]](#1), where trial design factors have been removed including AR1 adjustment of the residuals.
## Loblolly pine data set
The Loblolly pine data example was provided as a part of the genomic selection benchmark data collection in Genetics by Resende et al. [[2]](#2). The issue can be found [here](https://academic.oup.com/genetics/article/190/4/1503/6064084).

# Example
Analysing the Scots pine data. All traits are used here and any individual with missing data is removed for simplicity. Note that different softwares can be used (here R).
```R
require(dplyr)
require(ggplot2)
library(forcats)
library(brms) # for analysis of PCs
library(nadiv) # for NRM calculation

######### adjusted phenotypic records
dataF264_v2 <- 
  read.table(file = "../data/data_F264.txt",sep=",",header = T) ## remember to check your path!
######### pedigree
pedigreeF264 <- read.table(file = "../data/pedigree_F264.txt",sep=",",header=T)
######### remove individuals with missing records and scale the data
dataF264_v2 <- dataF264_v2 %>%
  na.omit() %>%
  mutate_each_(list(~scale(.) %>% as.vector),
                                  vars = 3:ncol(dataF264_v2))
```
Now, the data is ready for PCA!
```R
results <- prcomp(dataF264_v2[,3:ncol(dataF264_v2)])
results$Genotype_id <- dataF264_v2$Genotype_id
# Access the proportion of variance explained by each PC
var_explained <- results$sdev^2 / sum(results$sdev^2)
# Create a data frame for plotting
plot_data <- data.frame(
  PC = as.factor(paste0("PC",1:length(var_explained))),
  VarianceExplained = var_explained*100
)
# Reorder the factor levels of PC based on VarianceExplained
plot_data <- plot_data %>%
  mutate(PC = fct_reorder(PC, VarianceExplained, .desc = TRUE))

# Create a bar plot using ggplot2
ggplot(plot_data, aes(x = PC, y = VarianceExplained)) +
  geom_bar(stat = "identity", fill = "steelblue",color = "black") +
  labs(
    x = "",
    y = "Variance explained (%)",
    title = ""
  ) +
  theme_minimal() + 
  theme(axis.text = element_text(size=11, color = "black"), axis.title = element_text(size = 11),
        axis.text.x = element_text(angle = 45))
# loading plot
# PC1 vs PC2
PC <- as.data.frame(results$rotation)
# fix labels
x_label <- sprintf("PC1 (%.1f%%)", plot_data$VarianceExplained[1])
y_label <- sprintf("PC2 (%.1f%%)", plot_data$VarianceExplained[2])
ggplot(PC, aes(x = PC1, y = PC2)) +
  geom_point(size=2) +
   theme(axis.text = element_text(size = 11,color = "black"),legend.text = element_text(size = 11)) +
  theme_minimal() +
  geom_text(
    label = rownames(PC),
    nudge_x=0.03, nudge_y=0.03,
    check_overlap=T,
color = "black",
  size = 3)  +
  labs(color='Trait category') +
  xlab(x_label) +
  ylab(y_label)
```
The code will produce a histogram of the amount of phenotypic variance explained for each principal component (PC):
![Histogram](https://github.com/jonhar97/Reduced_phenotype_MME/blob/main/R/histogram_example.png)
A figure of the PC1 - PC2 plane of factor loadings highlights the trait relationships: 
![Ordination](https://github.com/jonhar97/Reduced_phenotype_MME/blob/main/R/ordination_pc1_pc2_example.png)
Lets continue to analyse the PCs as response variables. Note that I only use the first three PCs here for simplicity. Keep all if you dont want to loose any information. Furthermore, I use Rstan via BRMS package [[3]](#3) for the LMM analysis, but another software package can easily be used. Please see [[4]](#4) for additional examples (shown in Table 3).
```R
ped_tmp <- data.frame(Genotype_id = dataF264_v2$Genotype_id, sire = dataF264_v2$Dad_id, dam = dataF264_v2$Mum_id) %>% 
  prepPed()
Amat <- as.matrix(makeA(ped_tmp)) # Create A matrix

bf_PC1 <- bf(PC1 ~ 1 + (1 | a | gr(Genotype_id, cov = Amat)))
bf_PC2 <- bf(PC2 ~ 1 + (1 | a | gr(Genotype_id, cov = Amat)))
bf_PC3 <- bf(PC3 ~ 1 + (1 | a | gr(Genotype_id, cov = Amat)))

############## PC1
brms_m1.pc1 <- brm(
  bf_PC1, 
  data = dataF264_v2,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = 6, 
  cores = 6
)
plot(brms_m1.pc1)
summary(brms_m1.pc1)
mcmc_plot(brms_m1.pc1, type = "acf")
############## PC2
brms_m1.pc2 <- brm(
  bf_PC2, 
  data = dataF264_v2,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = 6, 
  cores = 6
)
plot(brms_m1.pc2)
summary(brms_m1.pc2)
mcmc_plot(brms_m1.pc2, type = "acf")
############## PC3
brms_m1.pc3 <- brm(
  bf_PC3, 
  data = dataF264_v2,
  data2 = list(Amat = Amat),
  family = gaussian(),
  chains = 6, 
  cores = 6
)
plot(brms_m1.pc3)
summary(brms_m1.pc3)
mcmc_plot(brms_m1.pc3, type = "acf")

###################### heritability ###############################
v_individual <- (VarCorr(brms_m1.pc1, summary = FALSE)$Genotype_id$sd)^2
v_r <- (VarCorr(brms_m1.pc1, summary = FALSE)$residual$sd)^2
h.pc1 <- as.mcmc(v_individual / (v_individual + v_r))
summary(h.pc1)
v_individual <- (VarCorr(brms_m1.pc2, summary = FALSE)$Genotype_id$sd)^2
v_r <- (VarCorr(brms_m1.pc2, summary = FALSE)$residual$sd)^2
h.pc2 <- as.mcmc(v_individual / (v_individual + v_r))
summary(h.pc2)
v_individual <- (VarCorr(brms_m1.pc3, summary = FALSE)$Genotype_id$sd)^2
v_r <- (VarCorr(brms_m1.pc3, summary = FALSE)$residual$sd)^2
h.pc3 <- as.mcmc(v_individual / (v_individual + v_r))
summary(h.pc3)

plot(h.pc1)
plot(h.pc2)
plot(h.pc3)
```
The code will (among other things) produce density plots of the heritability and a trace plot of the corresponding MCMC:
![Heritavility](https://github.com/jonhar97/Reduced_phenotype_MME/blob/main/R/h2_pc2.png)
Now to obtain EBV on the original scale, we can backtransform EBV on the PC scale:
```R
############################## Backtransform
### df with EBVs for each tree
posterior_samples1 <- posterior_samples(brms_m1.pc1)
posterior_samples2 <- posterior_samples(brms_m1.pc2)
posterior_samples3 <- posterior_samples(brms_m1.pc3)
# Calculate the posterior mean for the Genotype_id effect
G_pc <- data.frame(PC1 = colMeans(posterior_samples1[, grep("^r_Genotype_id\\[", colnames(posterior_samples1))]), PC2 = colMeans(posterior_samples2[, grep("^r_Genotype_id\\[", colnames(posterior_samples2))]), PC3 = colMeans(posterior_samples3[, grep("^r_Genotype_id\\[", colnames(posterior_samples3))])) %>%
  as.matrix()

rotation <- as.matrix(results$rotation[,1:3])
mu <- colMeans(dataF264_v2[,c("Adj_Dia_14","Adj_Dia_26","Adj_Dia_47","Adj_Fstam_47","Adj_Ftopp_47","Adj_Gant_14","Adj_Gdia_26","Adj_Gdiagr130_14","Adj_Gvin_26","Adj_Gvingr130_14","Adj_Hjd_10",     "Adj_Hjd_14","Adj_Hjd_26","Adj_Hjd_47","Adj_Kvar_47","Adj_Lev_10","Adj_Lev_14","Adj_Lev_26","Adj_Lev_47","Adj_Rak_14","Adj_Sdbst_14","Adj_Sprant_14","Adj_Sprant_26","Adj_Vit_10","Adj_Vit_14","Adj_Vit_26","Adj_Vit_47")])
EBV.pc <- as.data.frame(G_pc %*% t(rotation) + mu)

```
# References
<a id="1">[1]</a> Dutkowski, G. W., Costa E Silva, J., Gilmour, A. R., Wellendorf, H., & Aguiar, A. (2006). Spatial analysis enhances modelling of a wide variety of traits in forest genetic trials. Canadian Journal of Forest Research, 36(7), 1851–1870. 

<a id="2">[2]</a> Resende, J. F. R., Muñoz, P., Resende, M. D. V., Garrick, D. J., Fernando, R. L., Davis, J. M., Jokela, E. J., Martin, T. A., Peter, G. F., & Kirst, M. (2012). Accuracy of genomic selection methods in a standard data set of loblolly pine (Pinus taeda L.). Genetics, 190(4), 1503–1510. 

<a id="3">[3]</a> Bürkner, P.-C. (2017). brms: An R Package for Bayesian Multilevel Models Using Stan. Journal of Statistical Software, 80(1), 1–28. 

<a id="4">[4]</a> Ahlinder, J., Hall, D., Sountama, M., & Sillanpӓӓ, M. J. (2024). Principal component analysis revisited: fast multi-trait genetic evaluations with smooth convergence. G3 Genes|Genomes|Genetics, Accepted.