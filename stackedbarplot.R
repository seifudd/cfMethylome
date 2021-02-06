
# library
library(ggplot2)
library(tidyverse)
library(dplyr)
library(reshape2)
library(matrixStats)
library(hrbrthemes)
library(gtsummary)
library(qvalue)
 
# create a dataset
# specie <- c(rep("sorgho" , 3) , rep("poacee" , 3) , rep("banana" , 3) , rep("triticum" , 3) )
# condition <- rep(c("normal" , "stress" , "Nitrogen") , 4)
# value <- abs(rnorm(12 , 0 , 15))
# data <- data.frame(specie,condition,value)

decon_output = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/11_set_MKJ5201/02-methylseq-analysis-pipeline/04-tissue-of-origin-prediction/02-10X_in_atleast_2_samples/cfDNA_methatlas_betavalues_lung_CpGs10X_atleast_2samples_deconv_output.csv")
decon_output

decon_output_long <- melt(decon_output,
        # ID variables - all the variables to keep but not split apart on
    id.vars=c("Cell_Tissue_Type"),
        # The source columns
    measure.vars=c("S1","S2","S3","S4","S5","S6","S7","S8"),
        # Name of the destination column that will identify the original
        # column that the measurement came from
    variable.name="IDs",
    value.name="Cell_Tissue_Type_proportion"
)
decon_output_long

# Stacked + percent
ggplot(decon_output_long, aes(fill=Cell_Tissue_Type, y=Cell_Tissue_Type_proportion, x=IDs)) + 
    geom_bar(position="fill", stat="identity") +
 # set axis limits in coord_cartesian
 coord_cartesian(ylim = c(0.75, 1))

###################################################################################################################
#
#
# COVID-19 cfDNA MethylSeq deconvolution plot
#
#
###################################################################################################################

decon_output = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/16_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249/02_methylseq_analysis_pipeline/02_tissue_of_origin_prediction/01_coverage_5Xbysample/covid19_CpGs5X_deconv_output_mod.csv")
decon_output

col.names = colnames(decon_output)[2:dim(decon_output)[2]]

decon_output_long <- melt(decon_output,
        # ID variables - all the variables to keep but not split apart on
    id.vars=c("Cell_Tissue_Type"),
        # The source columns
    measure.vars = c("mcD-009-H2711M-1d","mcD-010-H2711M-8d","mcD-011-H2711M-11d","mcD-012-J2710G-1d","mcD-013-J2710G-13d","mcD-014-J2710G-16d","mcD-015-M18O1V-4d","mcD-016-M18O1V-16d","mcD-017-O0407LG-1d","mcD-018-O0407LG-9d","mcD-019-O0407LG-11d","mcD-020-P1201G-1d","mcD-021-P1201G-7d","mcD-022-P1201G-16d","mcD-023-T1909KA-1d","mcD-024-T1909KA-7d","BB295-Control","BB312-Control","BB380-Control"),
        # Name of the destination column that will identify the original
        # column that the measurement came from
    variable.name="IDs",
    value.name="Cell_Tissue_Type_proportion"
)
decon_output_long = as_tibble(decon_output_long)
decon_output_long

# measure.vars=c("mcD-009","mcD-010","mcD-011","mcD-012","mcD-013","mcD-014","mcD-015","mcD-016","mcD-017","mcD-018","mcD-019","mcD-020","mcD-021","mcD-022","mcD-023","mcD-024")
# decon_output_longsrt <- decon_output_long %>% group_by(IDs)
# decon_output_longsrt <- decon_output_longsrt %>% arrange(desc(Cell_Tissue_Type_proportion), .by_group = TRUE)
# decon_output_longsrt

# library(randomcoloR)
# n = length(levels(as.factor(decon_output_long$Cell_Tissue_Type)))
# col = distinctColorPalette(k=n, runTsne=T)

library("Polychrome")
data(palette36)
col = as.vector(palette36)
n = length(levels(as.factor(decon_output_long$Cell_Tissue_Type)))
col = sample(col, size=n)

# Stacked + percent
ggplot(decon_output_long, aes(x=IDs, y=Cell_Tissue_Type_proportion, fill=fct_reorder(Cell_Tissue_Type,Cell_Tissue_Type_proportion, .desc=T))) + 
	ggtitle("cfDNA_methatlas_covid19_CpGs5X_individual_samples") +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=col) +
    guides(fill = guide_legend(title="", ncol=1, label.theme = element_text(size = 7, angle = 0)
)) +
    theme(axis.text.x = element_text(size = 7, angle = 90), plot.title = element_text(size = 7))

 # Stacked + percent + zoomed
ggplot(decon_output_long, aes(x=IDs, y=Cell_Tissue_Type_proportion, fill=fct_reorder(Cell_Tissue_Type,Cell_Tissue_Type_proportion, .desc=T))) + 
	ggtitle("cfDNA_methatlas_covid19_CpGs5X_individual_samples_zoomed_in_between_0_0.30") +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=col) +
    guides(fill = guide_legend(title="", ncol=1, label.theme = element_text(size = 7, angle = 0)
)) +
    theme(axis.text.x = element_text(size = 7, angle = 90), plot.title = element_text(size = 7)) +
 # set axis limits in coord_cartesian
 coord_cartesian(ylim = c(0, 0.30))
dev.off()

###################################################################################################################
#
#
# COVID-19 cfDNA MethylSeq deconvolution plot (controls, batch1, batch2, batch3)
#
#
###################################################################################################################

setwd("/data/NHLBI_BCB/Sean_MethylSeq/17_MKJ5433/02_methylseq_analysis_pipeline/01-tissue_of_origin_prediction/01_coverage_5Xbysample/")

decon_output = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/17_MKJ5433/02_methylseq_analysis_pipeline/01-tissue_of_origin_prediction/01_coverage_5Xbysample/covid19_batch3_CpGs5X_deconv_output_mod.csv")
decon_output
col.names = colnames(decon_output)[2:dim(decon_output)[2]]
col.names

jhu_samples = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/17_MKJ5433/02_methylseq_analysis_pipeline/01-tissue_of_origin_prediction/01_coverage_5Xbysample/Transplant_patients_list_modFS.csv")
jhu_samples

decon_output = decon_output %>% select(c("Cell_Tissue_Type",jhu_samples$ID))
decon_output
col.names = colnames(decon_output)[2:dim(decon_output)[2]]
col.names

decon_output_long <- melt(decon_output,
        # ID variables - all the variables to keep but not split apart on
    id.vars=c("Cell_Tissue_Type"),
        # The source columns
    measure.vars = col.names,
        # Name of the destination column that will identify the original
        # column that the measurement came from
    variable.name="IDs",
    value.name="Cell_Tissue_Type_proportion"
)
decon_output_long = as_tibble(decon_output_long)
decon_output_long

# measure.vars=c("mcD-009","mcD-010","mcD-011","mcD-012","mcD-013","mcD-014","mcD-015","mcD-016","mcD-017","mcD-018","mcD-019","mcD-020","mcD-021","mcD-022","mcD-023","mcD-024")
# decon_output_longsrt <- decon_output_long %>% group_by(IDs)
# decon_output_longsrt <- decon_output_longsrt %>% arrange(desc(Cell_Tissue_Type_proportion), .by_group = TRUE)
# decon_output_longsrt

# library(randomcoloR)
# n = length(levels(as.factor(decon_output_long$Cell_Tissue_Type)))
# col = distinctColorPalette(k=n, runTsne=T)

library("Polychrome")
data(palette36)
col = as.vector(palette36)
n = length(levels(as.factor(decon_output_long$Cell_Tissue_Type)))
col = sample(col, size=n)

# Monocytes   
# Neutrophils 
# Erythrocyte_progenitors 
# Vascular_endothelial_cells  
# Adipocytes  
# Hepatocytes 
# Left_atrium 
# Kidney  
# Lung_cells  
# Pancreatic_beta_cells   
# Pancreatic_acinar_cells 
# Pancreatic_duct_cells   
# Colon_epithelial_cells  
# Bladder 
# Head_and_neck_larynx   
# Upper_GI    
# Uterus_cervix   
# Cortical_neurons
# Breast
# Prostate    
# Thyroid
# B-cells
# CD4T-cells
# NK-cells
# CD8T-cells

col = c("#1CFFCE",
    "#C075A6",
    "#B10DA1", 
    "#BDCDFF", 
    "#C4451C", 
    "#782AB6", 
    "#FBE426",
    "#683B79",
    "#3B00FB",
    "#B5EFB5",
    "#90AD1C",
    "#1C8356",
    "#E4E1E3",
    "#3283FE",
    "#B00068",
    "#16FF32",
    "#85660D",
    "#1C7F93",
    "#325A9B",
    "#F8A19F",
    "#DEA0FD",
    "#2ED9FF",
    "#AAF400",
    "#FC1CBF",
    "#D85FF7")

# Stacked + percent
ggplot(decon_output_long, aes(x=IDs, y=Cell_Tissue_Type_proportion, fill=fct_reorder(Cell_Tissue_Type,Cell_Tissue_Type_proportion, .desc=T))) + 
    ggtitle("cfDNA_methatlas_covid19_CpGs5X_individual_samples") +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=col) +
    guides(fill = guide_legend(title="", ncol=1, label.theme = element_text(size = 7, angle = 0)
)) +
    theme(axis.text.x = element_text(size = 3.5, angle = 90), plot.title = element_text(size = 7))

 # Stacked + percent + zoomed
ggplot(decon_output_long, aes(x=IDs, y=Cell_Tissue_Type_proportion, fill=fct_reorder(Cell_Tissue_Type,Cell_Tissue_Type_proportion, .desc=T))) + 
    ggtitle("cfDNA_methatlas_covid19_CpGs5X_individual_samples_zoomed_in_between_0_0.30") +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=col) +
    guides(fill = guide_legend(title="", ncol=1, label.theme = element_text(size = 7, angle = 0)
)) +
    theme(axis.text.x = element_text(size = 3.5, angle = 90), plot.title = element_text(size = 7)) +
 # set axis limits in coord_cartesian
 coord_cartesian(ylim = c(0, 0.30))
dev.off()

#########################################################################################################################
#
# Plot individual COVID-19 samples cell type proportions by time point
#
#
#########################################################################################################################

# setwd("/data/NHLBI_BCB/Sean_MethylSeq/20_plot_deconvolution_results_over_time_by_sample")
setwd("/data/NHLBI_BCB/Sean_MethylSeq/20_plot_deconvolution_results_over_time_by_sample/00-set2_plots")

# decon_output = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/20_plot_deconvolution_results_over_time_by_sample/00-set2_plots/longitudinal_plots.csv")
decon_output = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/20_plot_deconvolution_results_over_time_by_sample/00-set2_plots/longitudinal_plots_proportions.csv")
decon_output

# decon_output_long = decon_output %>% pivot_longer(!c(Patient_ID,PostECMO_Vent_Days),names_to = "Cell_Tissue_Type", values_to = "cfDNA_Concentration")
decon_output_long = decon_output %>% pivot_longer(!c(Patient_ID,PostECMO_Vent_Days),names_to = "Cell_Tissue_Type", values_to = "cfDNA_Proportions")
decon_output_long

library("Polychrome")
data(palette36)
col = as.vector(palette36)
n = length(levels(as.factor(decon_output_long$Cell_Tissue_Type)))
col = sample(col, size=n)
col

col = c("#1CFFCE",
    "#C075A6",
    "#B10DA1", 
    "#BDCDFF", 
    "#C4451C", 
    "#782AB6", 
    "#FBE426",
    "#683B79",
    "#3B00FB",
    "#B5EFB5",
    "#90AD1C",
    "#1C8356",
    "#E4E1E3",
    "#3283FE")


patientIDs = unique(decon_output$Patient_ID)
patientIDs

for (i in 1:length(patientIDs)) {
    print(patientIDs[i])
    decon_output_long = decon_output %>% pivot_longer(!c(Patient_ID,PostECMO_Vent_Days),names_to = "Cell_Tissue_Type", values_to = "cfDNA_Concentration")
    decon_output_long = filter(decon_output_long, Patient_ID == patientIDs[i])
    decon_output_long = filter(decon_output_long, cfDNA_Concentration>0)
    # Plot again
    plot = ggplot(decon_output_long, aes(x=PostECMO_Vent_Days, y=log10(cfDNA_Concentration), fill=Cell_Tissue_Type)) + 
        ggtitle(patientIDs[i]) +
        geom_area(alpha=0.70, outline.type="lower" , size=.2, colour="black") + 
        scale_fill_manual(values=col) +
        theme(legend.position="none")

    #save plots as .pdf
    ggsave(plot, file=paste(patientIDs[i], ".pdf", sep=''), scale=2)
}


for (i in 1:length(patientIDs)) {
    print(patientIDs[i])
#   decon_output = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/20_plot_deconvolution_results_over_time_by_sample/00-set2_plots/longitudinal_plots.csv")
#   decon_output_long = decon_output %>% pivot_longer(!c(Patient_ID,PostECMO_Vent_Days),names_to = "Cell_Tissue_Type", values_to = "cfDNA_Proportions")
    doutput = filter(decon_output_long, Patient_ID == patientIDs[i])
#   decon_output_long = filter(decon_output_long, cfDNA_Proportions>0)
    # Plot again

    pdf(file=paste(patientIDs[i],"_cfDNA_Proportions.pdf", sep=""), width=90/25.4, height=60/25.4, paper="special", bg="white", fonts="Helvetica", colormodel = "cmyk", pointsize=6, useDingbats=FALSE)

    myplot = ggplot(doutput, aes(x=PostECMO_Vent_Days, y=cfDNA_Proportions, fill=Cell_Tissue_Type)) + 
        ggtitle(patientIDs[i]) +
        geom_area(alpha=0.70, outline.type="lower" , size=.2, colour="black") + 
        scale_fill_manual(values=col) +
        scale_y_continuous(breaks = c(0,1)) +
        theme(legend.position="none") + 
        ylab("") +
        xlab("")
#       theme(panel.spacing = unit(0, 'pt'), strip.text = element_text(family='Helvetica', size=2), plot.margin=grid::unit(c(1,1,1,1), "pt")) +
#       facet_wrap(vars(patientIDs[i]), scales='free_x', ncol=1, strip.position = "right")

    print(myplot)
    #save plots as .pdf
#    ggsave(plot, file=paste(patientIDs[i], ".pdf", sep=''), scale=2)
    dev.off()
}


values=as.vector(pals::polychrome())[c(1:3,9,5:8,4,10,20,12:19,28)]



###################################################################################################################
#
#
# COVID-19 cfDNA MethylSeq deconvolution plot (batch4)
#
#
###################################################################################################################

library("here")
here()

decon_output = read_csv(file=here("02_methylseq_analysis_pipeline","02_tissue_of_origin_prediction","01_coverage_5Xbysample","covid19_batch4_CpGs5X_deconv_output.csv"))
decon_output

col.names = colnames(decon_output)[2:dim(decon_output)[2]]
col.names

jhu_samples = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/17_MKJ5433/02_methylseq_analysis_pipeline/01-tissue_of_origin_prediction/01_coverage_5Xbysample/Transplant_patients_list_modFS.csv")
jhu_samples

decon_output = decon_output %>% select(c("Cell_Tissue_Type",jhu_samples$ID))
decon_output
col.names = colnames(decon_output)[2:dim(decon_output)[2]]
col.names

decon_output_long <- melt(decon_output,
        # ID variables - all the variables to keep but not split apart on
    id.vars=c("Cell_Tissue_Type"),
        # The source columns
    measure.vars = col.names,
        # Name of the destination column that will identify the original
        # column that the measurement came from
    variable.name="IDs",
    value.name="Cell_Tissue_Type_proportion"
)
decon_output_long = as_tibble(decon_output_long)
decon_output_long

# measure.vars=c("mcD-009","mcD-010","mcD-011","mcD-012","mcD-013","mcD-014","mcD-015","mcD-016","mcD-017","mcD-018","mcD-019","mcD-020","mcD-021","mcD-022","mcD-023","mcD-024")
# decon_output_longsrt <- decon_output_long %>% group_by(IDs)
# decon_output_longsrt <- decon_output_longsrt %>% arrange(desc(Cell_Tissue_Type_proportion), .by_group = TRUE)
# decon_output_longsrt

# library(randomcoloR)
# n = length(levels(as.factor(decon_output_long$Cell_Tissue_Type)))
# col = distinctColorPalette(k=n, runTsne=T)

library("Polychrome")
data(palette36)
col = as.vector(palette36)
n = length(levels(as.factor(decon_output_long$Cell_Tissue_Type)))
col = sample(col, size=n)

col = c("#1CFFCE",
    "#C075A6",
    "#B10DA1", 
    "#BDCDFF", 
    "#C4451C", 
    "#782AB6", 
    "#FBE426",
    "#683B79",
    "#3B00FB",
    "#B5EFB5",
    "#90AD1C",
    "#1C8356",
    "#E4E1E3",
    "#3283FE",
    "#B00068",
    "#16FF32",
    "#85660D",
    "#1C7F93",
    "#325A9B",
    "#F8A19F",
    "#DEA0FD",
    "#2ED9FF",
    "#AAF400",
    "#FC1CBF",
    "#D85FF7")

# Stacked + percent
ggplot(decon_output_long, aes(x=IDs, y=Cell_Tissue_Type_proportion, fill=fct_reorder(Cell_Tissue_Type,Cell_Tissue_Type_proportion, .desc=T))) + 
    ggtitle("cfDNA_methatlas_covid19_CpGs5X_individual_samples") +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=col) +
    guides(fill = guide_legend(title="", ncol=1, label.theme = element_text(size = 7, angle = 0)
)) +
    theme(axis.text.x = element_text(size = 3.5, angle = 90), plot.title = element_text(size = 7))

 # Stacked + percent + zoomed
ggplot(decon_output_long, aes(x=IDs, y=Cell_Tissue_Type_proportion, fill=fct_reorder(Cell_Tissue_Type,Cell_Tissue_Type_proportion, .desc=T))) + 
    ggtitle("cfDNA_methatlas_covid19_CpGs5X_individual_samples_zoomed_in_between_0_0.30") +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=col) +
    guides(fill = guide_legend(title="", ncol=1, label.theme = element_text(size = 7, angle = 0)
)) +
    theme(axis.text.x = element_text(size = 3.5, angle = 90), plot.title = element_text(size = 7)) +
 # set axis limits in coord_cartesian
 coord_cartesian(ylim = c(0, 0.30))
dev.off()


###################################################################################################################
#
#
# COVID-19 cfDNA MethylSeq deconvolution plot (review)
#
#
###################################################################################################################

library("here")
here()

decon_output = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/covid19_deconvolution_allsamples_output_01202021_FS.txt")
decon_output

decon_output = decon_output %>% select(c("Cell_Tissue_Type","mcD_126","mcD_202","mcD_203"))
decon_output
colnames(decon_output) = c("Cell_Tissue_Type","covid19","healthy_control_1","healthy_control_2")
decon_output
col.names = colnames(decon_output)[2:dim(decon_output)[2]]
col.names

decon_output_long <- melt(decon_output,
        # ID variables - all the variables to keep but not split apart on
    id.vars=c("Cell_Tissue_Type"),
        # The source columns
    measure.vars = col.names,
        # Name of the destination column that will identify the original
        # column that the measurement came from
    variable.name="IDs",
    value.name="Cell_Tissue_Type_proportion"
)
decon_output_long = as_tibble(decon_output_long)
decon_output_long

# measure.vars=c("mcD-009","mcD-010","mcD-011","mcD-012","mcD-013","mcD-014","mcD-015","mcD-016","mcD-017","mcD-018","mcD-019","mcD-020","mcD-021","mcD-022","mcD-023","mcD-024")
# decon_output_longsrt <- decon_output_long %>% group_by(IDs)
# decon_output_longsrt <- decon_output_longsrt %>% arrange(desc(Cell_Tissue_Type_proportion), .by_group = TRUE)
# decon_output_longsrt

# library(randomcoloR)
# n = length(levels(as.factor(decon_output_long$Cell_Tissue_Type)))
# col = distinctColorPalette(k=n, runTsne=T)

library("Polychrome")
data(palette36)
col = as.vector(palette36)
n = length(levels(as.factor(decon_output_long$Cell_Tissue_Type)))
col = sample(col, size=n)

col = c("#1CFFCE",
    "#C075A6",
    "#B10DA1", 
    "#BDCDFF", 
    "#C4451C", 
    "#782AB6", 
    "#FBE426",
    "#683B79",
    "#3B00FB",
    "#B5EFB5",
    "#90AD1C",
    "#1C8356",
    "#E4E1E3",
    "#3283FE",
    "#B00068",
    "#16FF32",
    "#85660D",
    "#1C7F93",
    "#325A9B",
    "#F8A19F",
    "#DEA0FD",
    "#2ED9FF",
    "#AAF400",
    "#FC1CBF",
    "#D85FF7")

# Stacked + percent
ggplot(decon_output_long, aes(x=IDs, y=Cell_Tissue_Type_proportion, fill=fct_reorder(Cell_Tissue_Type,Cell_Tissue_Type_proportion, .desc=T))) + 
    ggtitle("cfDNA_methatlas_covid19_vs_healthy_control") +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=col) +
    guides(fill = guide_legend(title="", ncol=1, label.theme = element_text(size = 7, angle = 0)
)) +
    theme(axis.text.x = element_text(size = 7, angle = 90), plot.title = element_text(size = 7))

 # Stacked + percent + zoomed
ggplot(decon_output_long, aes(x=IDs, y=Cell_Tissue_Type_proportion, fill=fct_reorder(Cell_Tissue_Type,Cell_Tissue_Type_proportion, .desc=T))) + 
    ggtitle("cfDNA_methatlas_covid19_vs_healthy_control_zoomed_in_between_0_0.30") +
    geom_bar(position="fill", stat="identity") +
    scale_fill_manual(values=col) +
    guides(fill = guide_legend(title="", ncol=1, label.theme = element_text(size = 7, angle = 0)
)) +
    theme(axis.text.x = element_text(size = 7, angle = 90), plot.title = element_text(size = 7)) +
 # set axis limits in coord_cartesian
 coord_cartesian(ylim = c(0, 0.30))
dev.off()

###########################################################################################################################
#
#
#
###########################################################################################################################

decon_output = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/covid19_deconvolution_allsamples_output_01202021_FS.txt")
decon_output

decon_output = decon_output %>% select(c("Cell_Tissue_Type","S1","S2","S3","mcD_202","mcD_216","mcD_215","mcD_199","mcD_218","mcD_196","mcD_214","mcD_193","mcD_185","mcD_206","mcD_205","mcD_201","mcD_203","mcD_204","mcD_207","mcD_208"))
decon_output
decon_output
decon_output = decon_output %>% mutate(mean = rowMeans(across(where(is.numeric))))
decon_output

p<-ggplot(data=decon_output, aes(x=reorder(Cell_Tissue_Type,mean), y=mean)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  xlab("Cell_Tissue_Type") + ylab("Mean_of_Proportions_across_all_healthy_controls") +
  theme(axis.text.x = element_text(size = 7.0, angle = 90))
p
dev.off() 




# decon_output_long = filter(decon_output_long, Patient_ID == "G2111GD")

# ggplot(decon_output_long, aes(x=PostECMO_Vent_Days, y=Cell_Tissue_Type_proportion, fill=Cell_Tissue_Type)) + 
#        ggtitle("G2111GD") +
#        geom_area() + 
#        scale_fill_manual(values=col)

####################################################################################################################################################
####################################################################################################################################################
####################################################################################################################################################
+ geom_text(data=charts.data, aes(x = year, y = percentage, label = paste0(percentage,"%")), size=4)

fill <- c("#5F9EA0", "#E1B378")
p4 <- p4 + scale_fill_manual(values=fill)
p4


title.theme = element_text(size = 15, face = "italic", colour = "red", angle = 0)
panel.background = element_blank()


library(corrr)
reference_atlas = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/10-tissue_of_origin_methylation_project/meth_atlas_old/reference_atlas.csv")
data = as.data.frame(reference_atlas)
datadf = data[,2:dim(data)[2]]
rs = correlate(datadf, method="pearson", diagonal=1)


reference_atlas = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/10-tissue_of_origin_methylation_project/meth_atlas_old/atlas.csv")
data = as.data.frame(reference_atlas)
datadf = data[,2:dim(data)[2]]
rs = correlate(datadf, method="pearson", diagonal=1)



# correlate(d, method = "spearman", diagonal = 1)

rs %>% fashion()

rs %>% rplot(shape = 20, legend = TRUE, cex=0.75)
dev.off()


rs %>%
  rearrange(method = "MDS", absolute = FALSE) %>%
  shave() %>% 
  rplot(shape = 15, colors = c("red", "green"))


col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white", "cyan", "#007FFF", "blue", "#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue")) 
col4 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "#7FFF7F", "cyan", "#007FFF", "blue", "#00007F"))
whiteblack <- c("white", "black")

library(corrplot)
res = cor(datadf)
corrplot(res, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)

## text labels rotated 45 degrees
corrplot(res, type = "lower", tl.col = "black", tl.srt = 45, tl.cex = 0.50)
dev.off()

corrplot.mixed(res, lower.col = "black", number.cex = .5, tl.cex = 0.50, tl.srt = 45, tl.col = "black")
dev.off()

## using these color spectra
corrplot(res, order = "hclust", addrect = 2, col = col1(100))


###########################################################################################################################################
#
# QC covid19 data - DMC/DMR analysis
# 1. covid19 vs. healthy control
# --volcano plot (summarize # of CpGs, qvalue < 0.05 [hypo or hyper methylated], pvalue < 0.05 [hypo or hyper methylated])
# --pca (colored by batch + category)
# --heatmap (with CpGs, qvalue < 0.05)
#
###########################################################################################################################################

library(qvalue)

group1 = "covid19"
group2 = "healthy_control"
feature = "CpG"

setwd("/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/00_healthycontrols_vs_covid19/")

singleCpG_results = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/00_healthycontrols_vs_covid19/healthycontrols_vs_covid19_betacoefficients_pval_df.txt")
singleCpG_results
singleCpG_results = singleCpG_results %>% filter(!is.na(pvalue))
singleCpG_results
singleCpG_results = singleCpG_results %>% mutate(qvalue = qvalue(pvalue)$qvalues)
singleCpG_results

singleCpG_results = singleCpG_results %>%
  mutate(
    q_or_p_significant = case_when(
      qvalue < 0.05 ~ "qvaluesig",
      pvalue < 0.05 ~ "pvaluesig",
      pvalue < 0.05 & qvalue < 0.05 ~ "q_and_p_sig",
      TRUE          ~ "not_sig"
    )
  )
singleCpG_results

ggplot(singleCpG_results, aes(beta, -log10(pvalue))) + geom_point(aes(col=q_or_p_significant)) + scale_color_manual(values=c("gray","lightblue","darkblue")) +
coord_cartesian(xlim = c(-1.0, 1.0), expand = TRUE) + 
geom_hline(aes(yintercept = 1.3), color="black", linetype="dashed", size=0.5)
ggsave(file=paste("cfMethylome_pvalue",feature,"_",group1,"_",group2,".png", sep=""))
dev.off()

############## Summarize CpG results

library(gtsummary)

# make dataset with a few variables to summarize
# stats <- metadata %>% select(Category,Total_Sequences,Total_Sequences_post_trimming,Reads_uniquely_aligned,Mapping_efficiency,Duplication_rate,Bisulfite_conversion_rate)

stats <- singleCpG_results %>% select(beta, pvalue, qvalue, q_or_p_significant)

table2 <- 
  tbl_summary(
    stats,
    by = q_or_p_significant, # split table by group
    statistic = all_continuous() ~ "{mean} ({min} - {max})",
    missing = "no" # don't list missing data separately
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() %>%
  as_tibble()

table2
write.table(table2, file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/00_healthycontrols_vs_covid19/healthycontrols_vs_covid19_CpG_summary_table.txt")

# heights %>% 
#  group_by(sex) %>%
#  summarize(average = mean(height), standard_deviation = sd(height))

############## PCA

healthycontrols_covid19_betavalues = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/00_healthycontrols_vs_covid19/healthycontrols_vs_covid19_betavals_granges_df.csv")

singleCpG_results_q0.05 = singleCpG_results %>% filter(qvalue<0.05)
singleCpG_results_q0.05

healthycontrols_covid19_betavalues_cpgs_q0.05 = healthycontrols_covid19_betavalues %>% semi_join(singleCpG_results_q0.05, by = c("CpGID" = "CpGID_chr_bp"))
healthycontrols_covid19_betavalues_cpgs_q0.05

phenodata = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/cfCOVID.sample.metafile.12162020.FS.txt")
phenodata
phenodata = phenodata %>% filter(`1_Covid_vs_Healthy` != "NA")
phenodata

# a = healthycontrols_covid19_betavalues %>% select(-X1,-chromosome,-startpos, -endpos, -strand) %>% pivot_longer(cols = everything(), names_to = "SeqID", values_to = "methylation") %>% inner_join(phenodata, by = "SeqID")
# a

healthycontrols_covid19_betavalues_cpgs_q0.05_mod = healthycontrols_covid19_betavalues_cpgs_q0.05[,6:dim(healthycontrols_covid19_betavalues_cpgs_q0.05)[2]]
healthycontrols_covid19_betavalues_cpgs_q0.05_mod
healthycontrols_covid19_betavalues_cpgs_q0.05_df = t(healthycontrols_covid19_betavalues_cpgs_q0.05_mod)

library("missMDA")
healthycontrols_covid19_betavalues_cpgs_q0.05_imputed = imputePCA(healthycontrols_covid19_betavalues_cpgs_q0.05_df)
healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_df = healthycontrols_covid19_betavalues_cpgs_q0.05_imputed$completeObs
healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_df[1:10,1:10]

pca = prcomp(healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_df)

eigs <- pca$sdev^2
e1=eigs[1] / sum(eigs) * 100
e2=eigs[2] / sum(eigs) * 100
e3=eigs[3] / sum(eigs) * 100
e4=eigs[4] / sum(eigs) * 100
e5=eigs[5] / sum(eigs) * 100
spc1=sprintf('PC1 (%% %.1f)',e1)
spc2=sprintf('PC2 (%% %.1f)',e2)
spc3=sprintf('PC3 (%% %.1f)',e3)
spc4=sprintf('PC4 (%% %.1f)',e4)
spc5=sprintf('PC5 (%% %.1f)',e5)

pcatibble = as_tibble(pca$x, rownames="SeqID")
pcatibble
pcatibble_pheno = pcatibble %>% inner_join(phenodata, by="SeqID")
pcatibble_pheno

pdf(paste("CpGs_qvalue0.05_first_five_PCs",feature,"_",group1,"_",group2,".pdf", sep=""))
pairs(pca$x[,1:5], main=paste("scatter plot matrix of Principal Components_first_five_PCs_",feature,"_",group1,"_",group2, sep=""), pch=21, bg=c('#e41a1c','#377eb8')[as.factor(pcatibble_pheno$"1_Covid_vs_Healthy")], cex.main=0.75)
dev.off()

pdf(paste("CpGs_qvalue0.05_PC1_PC2",feature,"_",group1,"_",group2,".pdf", sep=""))
# Add extra space to right of plot area; change clipping to figure
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(pca$x[,1], pca$x[,2], main=paste("scatter plot matrix of Principal Components_PC1.vs.PC2_",feature,"_",group1,"_",group2, sep=""), pch=16, col=c('#e41a1c','#377eb8')[as.factor(pcatibble_pheno$"1_Covid_vs_Healthy")], cex=1.20, cex.main=0.80, xlab=spc1, ylab=spc2)
# text(pca$x[,1], pca$x[,2], labels=paste(rownames(phenodata),phenodata$dx,sep="_"), pos=3, cex=0.30)
legend("topright", inset=c(-0.23,0), legend=c(group1,group2), pch=c(16,16), col=rev(levels(as.factor(c('#e41a1c','#377eb8')[as.factor(pcatibble_pheno$"1_Covid_vs_Healthy")]))), title="Group", cex=0.60)
dev.off()

############## Heatmap for specific group comparisons

col_fun = colorRamp2(seq(min(healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_df, na.rm=TRUE), max(healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_df, na.rm=TRUE), length = 3), c("green", "black", "red"), space = "LUV")

healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_tibble = as_tibble(healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_df, rownames="SeqID")
healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_tibble
healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_tibble_pheno = healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_tibble %>% inner_join(phenodata, by="SeqID")
healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_tibble_pheno

groups = factor(healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_tibble_pheno$"1_Covid_vs_Healthy", levels = c("HC", "Covid"))
groups
column_ha = HeatmapAnnotation(group = groups, 
    annotation_legend_param = list(title=c("group"), title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6), group=list(at=c("HC", "Covid"))), 
    col = list(group = c("HC" = "blue", "Covid" = "yellow")), 
    annotation_label = "group",
    annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

b = as.matrix(t(healthycontrols_covid19_betavalues_cpgs_q0.05_imputed_df))

ht = Heatmap(b, 
    col = col_fun,
    top_annotation = column_ha, 
    column_title = "Patient_Category", 
    row_title = "CpGs_Methylation", 
    heatmap_legend_param = list(title = "Methylation", title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6)), 
    column_title_gp = gpar(fontsize = 7), 
    row_title_gp = gpar(fontsize = 7),
    row_names_gp = gpar(fontsize = 8), 
    show_column_names = TRUE,
    show_row_names = FALSE,
    column_names_gp = gpar(fontsize = 3),  
    cluster_columns = TRUE, 
    cluster_rows = TRUE,
    width = ncol(b)*unit(1.8, "mm"), 
    height = ncol(b)*unit(1.5, "mm"))

calc_ht_size = function(ht, unit = "inch") {
    pdf(NULL)
    ht = draw(ht)
    w = ComplexHeatmap:::width(ht)
    w = convertX(w, unit, valueOnly = TRUE)
    h = ComplexHeatmap:::height(ht)
    h = convertY(h, unit, valueOnly = TRUE)
    dev.off()

    c(w, h)
}

size = calc_ht_size(ht)
size
pdf("Rplots.pdf", width=size[1], height=size[2])
ht
dev.off()

############## Distribution of methylation across samples

ggplot(aes(x=methylation, colour=Library_ID), data=a) + 
    geom_density(show.legend=FALSE)
dev.off()

a
ggplot(aes(y=methylation, colour=Library_ID), data=a) + 
    geom_boxplot(aes(colour=Library_ID, y=methylation), show.legend=FALSE)
dev.off()


############## Summarize DMR results

regionp = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/00_healthycontrols_vs_covid19/01-combp/cpv_genomewide_healthycontrols_vs_covid19.regions-p.bed")
regionp

regionpt = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/00_healthycontrols_vs_covid19/01-combp/cpv_genomewide_healthycontrols_vs_covid19.regions-t.bed")
regionpt

table2 <- 
  tbl_summary(
    regionpt,
    statistic = all_continuous() ~ "{mean} ({min} - {max})",
    missing = "no" # don't list missing data separately
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() %>%
  as_tibble()

table2

regionpt = regionpt %>% mutate(hypo_hyper = case_when(`t.sum` < 0.0 ~ 1, `t.sum` > 0.0 ~ 0, TRUE ~ NA_real_))

regionpt_annotated = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/00_healthycontrols_vs_covid19/01-combp/00-DMRs_annotated/cpv_genomewide_healthycontrols_vs_covid19.regions-t.sorted.annotated.bed")
regionpt_annotated


########################################################################################################################################
#
#
# plot deconvolution results (sum across all samples by Cell_Tissue_Type)
# purpose:  which Cell_Tissue_Type to include as covariates in the linear model? 
# 
#
########################################################################################################################################

props = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/covid19_deconvolution_allsamples_output_11302020_FS.txt")
sumprops = props %>% summarise(across(c("Monocytes_EPIC":"Uterus_cervix"), ~sum(.x, na.rm = TRUE)))
sumprops = t(sumprops)
sumprops = as_tibble(as.data.frame(sumprops)) %>% add_column(rownames(sumprops), .before="V1") %>% arrange(V1)
colnames(sumprops) = c("Cell_Tissue_Type", "sum_proportions_all_samples")
sumprops

p<-ggplot(data=sumprops, aes(x=reorder(Cell_Tissue_Type,sum_proportions_all_samples), y=sum_proportions_all_samples)) +
  geom_bar(stat="identity", show.legend=FALSE) +
  geom_hline(aes(yintercept = 4.0), color="black", linetype="dashed",  size=0.5) +
  xlab("Cell_Tissue_Type") + ylab("Sum_of_Proportions_across_all_Samples") +
  theme(axis.text.x = element_text(size = 7.0, angle = 90))
p
dev.off() 

########################################################################################################################################
#
# all samples, all CpGs
# 1. heat map of CpGs methylation across all samples (by group)
# 2. circos plot of CpGs methylation across all samples (by group)
#
########################################################################################################################################

coverage_all_samples = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/00_methylation_estimates_and_coverage_by_chromosome_across_all_samples_all_CpGs/all_samples_covid19_hc_influenza_coverage_granges_df.csv")
coverage_all_samples

samples = colnames(coverage_all_samples)[6:dim(coverage_all_samples)[2]]
length(samples)

coverage_all_samples_v2 = coverage_all_samples %>% mutate(covfilter=rowSums(.[samples] >= 5) >= 112)
coverage_all_samples_v2
coverage_all_samples_v3 = coverage_all_samples_v2 %>% filter(covfilter=="TRUE")
coverage_all_samples_v3
keepLoci = as_tibble(coverage_all_samples_v3$CpGID)
colnames(keepLoci) = "CpGID"
keepLoci

betavalues_all_samples = read_csv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/00_methylation_estimates_and_coverage_by_chromosome_across_all_samples_all_CpGs/all_samples_covid19_hc_influenza_betavals_granges_df.csv")
betavalues_all_samples

betavalues_all_samples_KeepLoci = betavalues_all_samples %>% semi_join(keepLoci, by = "CpGID")
betavalues_all_samples_KeepLoci


hc = phenodata %>% filter(Category=="HC")
hc = hc$SeqID
length(hc)

icu = phenodata %>% filter(Category=="ICU")
icu = icu$SeqID
length(icu)

nonicu = phenodata %>% filter(Category=="Non-ICU")
nonicu = nonicu$SeqID
length(nonicu)

influenza_rsv = phenodata %>% filter(Category=="Influenza & RSV")
influenza_rsv = influenza_rsv$SeqID
length(influenza_rsv)

mild_covid19 = phenodata %>% filter(Category=="Mild COVID-19")
mild_covid19 = mild_covid19$SeqID
length(mild_covid19)

transplant = phenodata %>% filter(Category=="COVID-19 + Transplant")
transplant  = transplant$SeqID
length(transplant)


betavalues_all_samples_KeepLoci_v2 = betavalues_all_samples_KeepLoci %>% 
  mutate(beta_mean_hc=rowMeans(.[hc], na.rm=TRUE), beta_mean_icu=rowMeans(.[icu], na.rm=TRUE), beta_mean_nonicu=rowMeans(.[nonicu], na.rm=TRUE), beta_mean_influenza_rsv=rowMeans(.[influenza_rsv], na.rm=TRUE), beta_mean_mild_covid19=rowMeans(.[mild_covid19], na.rm=TRUE), beta_mean_transplant=rowMeans(.[transplant], na.rm=TRUE))
betavalues_all_samples_KeepLoci_v2

groupbeta_colnames = c("beta_mean_hc", "beta_mean_icu", "beta_mean_nonicu", "beta_mean_influenza_rsv", "beta_mean_mild_covid19", "beta_mean_transplant")
betavalues_all_samples_KeepLoci_v3 = betavalues_all_samples_KeepLoci_v2 %>% mutate(beta_var=rowVars(as.matrix(.[groupbeta_colnames])))
betavalues_all_samples_KeepLoci_v3

summary(betavalues_all_samples_KeepLoci_v3$beta_var)

ggplot(betavalues_all_samples_KeepLoci_v3, aes(x=beta_var)) + 
 geom_histogram(aes(y=..density..), colour="black", fill="white")+
 geom_density(alpha=.2, fill="#FF6666") +
 xlab("CpG_methylation_variance_across_mean_methylation_by_category")
dev.off()

betavalues_all_samples_filtered_by_variance = betavalues_all_samples_KeepLoci_v3 %>% filter(beta_var>=0.06) # beta_var>=0.03, 2,573 CpGs
betavalues_all_samples_filtered_by_variance
betavalues_all_samples_filtered_by_variance = betavalues_all_samples_filtered_by_variance %>% select(-c("beta_mean_hc", "beta_mean_icu", "beta_mean_nonicu", "beta_mean_influenza_rsv", "beta_mean_mild_covid19", "beta_mean_transplant", "beta_var"))
betavalues_all_samples_filtered_by_variance_df = betavalues_all_samples_filtered_by_variance
# betavalues_all_samples_filtered_by_variance_df = betavalues_all_samples_filtered_by_variance[,6:dim(betavalues_all_samples_filtered_by_variance)[2]]
betavalues_all_samples_filtered_by_variance_df
# betavalues_all_samples_filtered_by_variance_df = betavalues_all_samples_filtered_by_variance_df %>% select(-c("mcD_92","mcD_93"))
betavalues_all_samples_filtered_by_variance_df

samples = as_tibble(colnames(betavalues_all_samples_filtered_by_variance_df))
colnames(samples) = "SeqID"
samples
phenodata_mod = samples %>% inner_join(phenodata, by = "SeqID")
phenodata_mod
phenodata_mod = phenodata_mod %>% filter(Category != "No information")
phenodata_mod

phenodata_subset = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/cfCOVID.sample.metafile.v2.12292020.FS.txt")
phenodata_subset

col_fun = colorRamp2(seq(min(betavalues_all_samples_filtered_by_variance_df, na.rm=TRUE), max(betavalues_all_samples_filtered_by_variance_df, na.rm=TRUE), length = 3), c("green", "black", "red"), space = "LUV")

groups = factor(phenodata_mod$Category, levels = c("HC", "Influenza & RSV", "Mild COVID-19", "ICU", "Non-ICU", "COVID-19 + Transplant"))
groups
column_ha = HeatmapAnnotation(group = groups, 
    annotation_legend_param = list(title=c("group"), title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6), group=list(at=c("HC", "Influenza & RSV", "Mild COVID-19", "Non-ICU", "ICU", "COVID-19 + Transplant"))), 
    col = list(group = c("COVID-19 + Transplant" = "brown", "Influenza & RSV" = "orange", "HC" = "blue", "Non-ICU" = "yellow", "Mild COVID-19" = "pink", "ICU" = "purple")), 
    annotation_label = "group",
    annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

b = as.matrix(betavalues_all_samples_filtered_by_variance_df)

ht = Heatmap(b, 
    col = col_fun,
    top_annotation = column_ha, 
    column_title = "Patient_Category", 
    row_title = "CpGs_Methylation", 
    heatmap_legend_param = list(title = "Methylation", title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6)), 
    column_title_gp = gpar(fontsize = 7), 
    row_title_gp = gpar(fontsize = 7),
    row_names_gp = gpar(fontsize = 8), 
    show_column_names = TRUE,
    show_row_names = FALSE,
    column_names_gp = gpar(fontsize = 3),  
    cluster_columns = TRUE, 
    width = ncol(b)*unit(0.8, "mm"), 
    height = ncol(b)*unit(0.5, "mm"))

calc_ht_size = function(ht, unit = "inch") {
    pdf(NULL)
    ht = draw(ht)
    w = ComplexHeatmap:::width(ht)
    w = convertX(w, unit, valueOnly = TRUE)
    h = ComplexHeatmap:::height(ht)
    h = convertY(h, unit, valueOnly = TRUE)
    dev.off()

    c(w, h)
}

size = calc_ht_size(ht)
size
pdf("Rplots.pdf", width=size[1], height=size[2])
ht
dev.off()

samples_cluster_order = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/00_methylation_estimates_and_coverage_by_chromosome_across_all_samples_all_CpGs/order_of_samples_heatmap.txt")
samples_cluster_order

samples_cluster_order = samples_cluster_order %>% inner_join(phenodata, by = "SeqID")
samples_cluster_order

########################################################################################################################################
#
# Create summary table of metadata/phenodata
#
########################################################################################################################################

library(gtsummary)

metadata = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/cfCOVID.sample.metafile.12162020.FS.txt")
metadata

# make dataset with a few variables to summarize
# stats <- metadata %>% select(Category,Total_Sequences,Total_Sequences_post_trimming,Reads_uniquely_aligned,Mapping_efficiency,Duplication_rate,Bisulfite_conversion_rate)
stats <- metadata %>% select(Category,Status,Total_Sequences,Total_Sequences_post_trimming,Reads_uniquely_aligned,Mapping_efficiency,Duplication_rate)

table2 <- 
  tbl_summary(
    stats,
    by = Category, # split table by group
    statistic = all_continuous() ~ "{mean} ({min} - {max})",
    missing = "no" # don't list missing data separately
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() %>%
  as_tibble()

table2
write.csv(table2, file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/stats.csv")

stats <- metadata %>% select("1_Covid_vs_Healthy", "2_ICU_vs_non_ICU", "3_recovered_vs_deceased", "5_Discharge_Date_off_ECMO")

table2 <- 
  tbl_summary(
    stats,
    statistic = all_categorical() ~ "{n} ({p}%)",
    missing = "no" # don't list missing data separately
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() %>%
  as_tibble()

table2
write.csv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/stats_breakdown_by_category_for_comparisons.csv", table2)

table2 <- 
  tbl_summary(
    stats,
    by = "1_Covid_vs_Healthy", # split table by group
    statistic = all_categorical() ~ "{n} ({p}%)",
    missing = "no" # don't list missing data separately
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() %>%
  as_tibble()

table2
write.csv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/stats_breakdown_by_covid19_vs_healthycontrol_for_comparisons.csv", table2)


stats <- metadata %>% select("Category", "Status", "5_Discharge_Date_off_ECMO")

table2 <- 
  tbl_summary(
    stats,
    by = "5_Discharge_Date_off_ECMO", # split table by group
    statistic = all_categorical() ~ "{n} ({p}%)",
    missing = "no" # don't list missing data separately
  ) %>%
  add_n() %>% # add column with total number of non-missing observations
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() %>%
  as_tibble()

table2
write.csv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/stats_breakdown_by_category_for_comparisons_disease_severity.csv", table2)


########################################################################################################################################
#
#
#
########################################################################################################################################

#########################################################################################################################
#
#
# Heatmap of deconvolution data of all samples using absolute concentrations or proportions
#
#
#########################################################################################################################

library(ComplexHeatmap)
library(pheatmap)
library(circlize)

setwd("/data/NHLBI_BCB/Sean_MethylSeq/28_heatmap_deconvolution_results")

# metadata = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/cfCOVID.sample.metafile.v2.12292020.FS.txt")
# metadata

# metadata = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/cfCOVID.sample.metafile.12162020.FS.txt")
# metadata

# metadata = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/28_heatmap_deconvolution_results/Transplant_COVID-19.txt")
# metadata

metadata = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/24_MKJ5564_MKJ5433_MKJ5389_MKJ5390_MKJ5391_MKJ5392_MKJ5249_MKJ5201/cfCOVID.sample.metafile.v3.01112021.FS.txt")
metadata

influenza_and_rsv = metadata %>% filter(Category_broad=="Influenza & RSV")
influenza_and_rsv
covid19_and_hc = metadata %>% filter(`1_Covid_vs_Healthy`!="NA")
covid19_and_hc

metadata = rbind(influenza_and_rsv,covid19_and_hc)
metadata

props = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/28_heatmap_deconvolution_results/Absolute_quantification_of_allsamples_output_12152020_FS.txt")
props
# props = read_tsv(file="/data/NHLBI_BCB/Sean_MethylSeq/28_heatmap_deconvolution_results/deconvolution_proportions_of_allsamples_output_12162020_FS.txt")
# batch1 = metadata %>% filter(Batch==1)
# props = props %>% filter(!(Library_ID %in% batch1$SeqID))

props = props %>% inner_join(metadata, by = "SeqID")
props

propsdf = props %>% select(Monocytes_EPIC:Uterus_cervix)
propsdf

propsdf = data.frame(t(propsdf))
propsdf[propsdf == 0] <- 1
propsdf = log10(propsdf)
colnames(propsdf) = props$SeqID
propsdf = as.matrix(propsdf)
propsdf[1:10,1:10]

# col_fun = colorRamp2(seq(min(propsdf), max(propsdf), length = 3), c("green", "black", "red"), space = "LUV")

col_fun = colorRamp2(seq(min(propsdf), max(propsdf), length = 3), c("blue", "#EEEEEE", "red"), space = "RGB")

################################################################################################################################################

groups = factor(props$Category_broad, levels = c("HC", "Influenza & RSV", "COVID-19"))
groups
column_ha = HeatmapAnnotation(group = groups, 
   annotation_legend_param = list(title=c("group"), title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 15), group=list(at=c("HC", "Influenza & RSV", "COVID-19"))), 
    col = list(group = c("Influenza & RSV" = "orange", "HC" = "blue", "COVID-19" = "black")), 
    annotation_label = "group",
    annotation_name_gp = gpar(fontsize = 15, fontface = "bold"))

################################################################################################################################################

################################################################################################################################################

# props$Category[props$Category == "COVID-19 + Transplant"] = "ICU"
# groups = factor(props$Category, levels = c("HC", "Influenza & RSV", "Mild COVID-19", "ICU", "Non-ICU"))
# groups
# column_ha = HeatmapAnnotation(group = groups, 
#    annotation_legend_param = list(title=c("group"), title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6), group=list(at=c("HC", "Influenza & RSV", "Mild COVID-19", "Non-ICU", "ICU"))), 
#    col = list(group = c("Influenza & RSV" = "orange", "HC" = "blue", "Non-ICU" = "yellow", "Mild COVID-19" = "pink", "ICU" = "purple")), 
#    annotation_label = "group",
#    annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

################################################################################################################################################

# groups = factor(props$Category, levels = c("HC", "Influenza & RSV", "Mild COVID-19", "ICU", "Non-ICU", "COVID-19 + Transplant"))
# groups
# column_ha = HeatmapAnnotation(group = groups, 
#    annotation_legend_param = list(title=c("group"), title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6), group=list(at=c("HC", "Influenza & RSV", "Mild COVID-19", "Non-ICU", "ICU", "COVID-19 + Transplant"))), 
#    col = list(group = c("COVID-19 + Transplant" = "brown", "Influenza & RSV" = "orange", "HC" = "blue", "Non-ICU" = "yellow", "Mild COVID-19" = "pink", "ICU" = "purple")), 
#    annotation_label = "group",
#    annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))

################################################################################################################################################

################################################################################################################################################
# groups = factor(props$Category, levels = c("Healthy Control", "Non-ICU COVID-19", "COVID-19 + Transplant"))
# groups
#column_ha = HeatmapAnnotation(group = groups, 
#   annotation_legend_param = list(title=c("group"), title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6), group=list(at=c("Healthy Control", "Non-ICU COVID-19", "COVID-19 + Transplant"))), 
#    col = list(group = c("COVID-19 + Transplant" = "brown", "Healthy Control" = "blue", "Non-ICU COVID-19" = "yellow")), 
#    annotation_label = "group",
#    annotation_name_gp = gpar(fontsize = 7, fontface = "bold"))
################################################################################################################################################

### small number of columns/samples and/or features/rows
ht = Heatmap(propsdf, 
    col = col_fun,
    top_annotation = column_ha, 
    column_title = "cfDNA_deconvolution_absolute_concentrations", 
    row_title = "Cell_Tissue_Type", 
    heatmap_legend_param = list(title = "log10_cfDNA_Concentration", title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6)), 
    column_title_gp = gpar(fontsize = 7), 
    row_title_gp = gpar(fontsize = 7),
    row_names_gp = gpar(fontsize = 6), 
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 5),  
    cluster_columns = TRUE, 
    width = ncol(propsdf)*unit(9, "mm"), 
    height = ncol(propsdf)*unit(10, "mm"))

calc_ht_size = function(ht, unit = "inch") {
    pdf(NULL)
    ht = draw(ht)
    w = ComplexHeatmap:::width(ht)
    w = convertX(w, unit, valueOnly = TRUE)
    h = ComplexHeatmap:::height(ht)
    h = convertY(h, unit, valueOnly = TRUE)
    dev.off()

    c(w, h)
}

size = calc_ht_size(ht)
size
pdf("Rplots.pdf", width=size[1], height=size[2])
ht
dev.off()

######################

ht = Heatmap(propsdf, 
    col = col_fun,
    top_annotation = column_ha, 
    column_title = "cfDNA_deconvolution_absolute_concentrations", 
    row_title = "Cell_Tissue_Type", 
    heatmap_legend_param = list(title = "log10_cfDNA_Concentration", title_gp = gpar(fontsize = 6, fontface = "bold"), labels_gp = gpar(fontsize = 6)), 
    column_title_gp = gpar(fontsize = 7), 
    row_title_gp = gpar(fontsize = 7),
    row_names_gp = gpar(fontsize = 6), 
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 3),  
    cluster_columns = TRUE, 
    width = ncol(propsdf)*unit(0.8, "mm"), 
    height = ncol(propsdf)*unit(0.5, "mm"))

calc_ht_size = function(ht, unit = "inch") {
    pdf(NULL)
    ht = draw(ht)
    w = ComplexHeatmap:::width(ht)
    w = convertX(w, unit, valueOnly = TRUE)
    h = ComplexHeatmap:::height(ht)
    h = convertY(h, unit, valueOnly = TRUE)
    dev.off()

    c(w, h)
}

size = calc_ht_size(ht)
size
pdf("Rplots.pdf", width=size[1], height=size[2])
ht
dev.off()

######################

ht = Heatmap(propsdf, 
    col = col_fun,
    top_annotation = column_ha, 
    column_title = "cfDNA_deconvolution_absolute_concentrations", 
    row_title = "Cell_Tissue_Type", 
    heatmap_legend_param = list(title = "log10_cfDNA_Concentration", title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13)), 
    column_title_gp = gpar(fontsize = 13, fontface = "bold"), 
    row_title_gp = gpar(fontsize = 13, fontface = "bold"),
    row_names_gp = gpar(fontsize = 15, fontface = "bold"), 
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 8, fontface = "bold"),  
    cluster_columns = TRUE, 
    width = ncol(propsdf)*unit(4, "mm"), 
    height = ncol(propsdf)*unit(3, "mm"))

calc_ht_size = function(ht, unit = "inch") {
    pdf(NULL)
    ht = draw(ht)
    w = ComplexHeatmap:::width(ht)
    w = convertX(w, unit, valueOnly = TRUE)
    h = ComplexHeatmap:::height(ht)
    h = convertY(h, unit, valueOnly = TRUE)
    dev.off()

    c(w, h)
}

size = calc_ht_size(ht)
size
pdf("Rplots.pdf", width=size[1], height=size[2])
ht
dev.off()


########################################################################################################################################
#
#
#
########################################################################################################################################

# metadata <- data.frame(props$Disease_type,row.names=colnames(propsdf))
# colnames(metadata) <- c("group")

# paletteLength = 150
# col = colorRampPalette( c("green", "black", "red"), space="rgb")(paletteLength)

pheatmap(propsdf,
    filename = paste("pheatmap","_","cfDNA_deconvolution_absolute_concentrations",".pdf",sep=""), 
    fontsize_row = 5, 
    fontsize_col = 3, 
    fontsize = 3, 
    show_rownames = T, 
    scale="none", 
    cluster_rows=T, 
    cluster_cols=T, 
    col=col, 
    cellwidth=2, 
    cellheight=15, 
    main='cfDNA_deconvolution', 
    border_color=NA, 
    annotation_col = metadata)

# annotation_row = data.frame(gene_direction_ibrutinib = as.vector(statres_gsigs$Direction.of.change))
# ann_colors = list(gene_direction_ibrutinib=c(UP="yellow", DN="blue"), timepoint=c(pre="orange", sixm="purple"))
# rownames(annotation_row) = make.names(statres_gsigs$genesymbol, unique=T)
# annotation_row=annotation_row, annotation_colors=ann_colors

# x=reorder(modules,corr)

chr22_betavals_granges_df = chr22_betavals_granges_df %>% mutate(Dx = case_when(Dx == "Covid" ~ 1, Dx == "HC" ~ 0, TRUE ~ NA_real_))

chr22_betavals_granges_df %>% select(c("Library_ID",starts_with("chr"),"Monocytes_EPIC","Neutrophils_EPIC","Erythrocyte_progenitors","Erythrocyte_progenitors","Adipocytes","Hepatocytes","Vascular_endothelial_cells","Left_atrium","1_Covid_vs_Healthy"))

chr22_betavals_granges_df = chr22_betavals_granges_df %>% select(c("Library_ID",starts_with("chr"),"Monocytes_EPIC","Neutrophils_EPIC","Erythrocyte_progenitors","Erythrocyte_progenitors","Adipocytes","Hepatocytes","Vascular_endothelial_cells","Left_atrium","1_Covid_vs_Healthy"))

chr22_betavals_granges_df = chr22_betavals_granges_df %>% mutate(Dx = case_when("1_Covid_vs_Healthy" == "Covid" ~ 1, "1_Covid_vs_Healthy" == "HC" ~ 0, TRUE ~ NA_real_))
