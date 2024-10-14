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
![Histogram of variance explained](https://github.com/jonhar97/Reduced_phenotype_MME/blob/main/R/histogram_example.png)


# References
<a id="1">[1]</a> Dutkowski, G. W., Costa E Silva, J., Gilmour, A. R., Wellendorf, H., & Aguiar, A. (2006). Spatial analysis enhances modelling of a wide variety of traits in forest genetic trials. Canadian Journal of Forest Research, 36(7), 1851–1870. 

<a id="2">[2]</a> Resende, J. F. R., Muñoz, P., Resende, M. D. V., Garrick, D. J., Fernando, R. L., Davis, J. M., Jokela, E. J., Martin, T. A., Peter, G. F., & Kirst, M. (2012). Accuracy of genomic selection methods in a standard data set of loblolly pine (Pinus taeda L.). Genetics, 190(4), 1503–1510. 