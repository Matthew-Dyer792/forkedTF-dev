## ForkedTF R_project 
# Denis Thieffry, May 15, 2023
# Project implemented on a Mac laptop  M1 (Arm64)
# With MacOS 12.6.2 (Monterey)
# With # R 4.3.0 - RSudio 2023.06.0+421 
# https://github.com/benoukraflab/forkedTF
# https://github.com/benoukraflab/forkedTF/wiki/ForkedTF
# Packages required to run forkedTF
install.packages("devtools")
library("devtools")
devtools::install_github("https://github.com/benoukraflab/forkedTF")
devtools::install_github("Matthew-Dyer792/tfregulomer-dev")
install.packages("BiocManager")
library("BiocManager")
BiocManager::install("ggplot2")
BiocManager::install("gridExtra")
BiocManager::install("gridGraphics")
BiocManager::install("ggseqlogo")
BiocManager::install("cowplot")
BiocManager::install("cowplot")
BiocManager::install("GenomicRanges")
#
# Set the default pathway for this project
main_dir <- "/Users/denis/Documents/R/RStudio_Projects/ForkedTF"
setwd(main_dir)
db_path <- file.path(main_dir, "TFregulomeR_database_2.1/tfregulome.sqlite")
# Load the TFregulomeR library and set the working directory for further analyses
# library("TFregulomeR")
# TFregulome_url = "https://methmotif.org/"
# all_records <- dataBrowser()
# head(all_records)
#
## ForkedTF Case study 1: JUND in HepG2 cells
# Create a subdirectory for storing the output of the HepG2 JUND analysis
sub_dir <- "HepG2_JUND"
dir.create(file.path(main_dir, sub_dir))
setwd(file.path(main_dir, sub_dir))
## MiniCofactor Report
# This function helps in the exploration of binding partners in a cell line. Input the mainTF, cell line of interest and the cobinding_threshold to generate a PDF report of the binding partners.
library(forkedTF)
## 
miniCofactorReport( TF = "JUND", cell = "HepG2" ) # This takes several minutes...
# In addition to finding the factor with the highest peak overlap, using the parameter filterBy="q.significance" we can compute a -log10(Adjusted P-value) from an enrichment test as implemented in https://github.com/remap-cisreg/ReMapEnrich
miniCofactorReport(TF = "JUND",cell = "HepG2", filterBy="q.significance") # This takes several minutes...
#
## FPWM creation and plot
# Use the createFPWM function to extract the motif, from empirical datasets, that a TF uses when binding with a partner TF. plotFPWM helps in visualizing the FPWM.
fpwm <- createFPWM(mainTF ="JUND",
                   partners = c("ATF7","ATF2","JUN","FOSL2","FOS"),
                   cell = "HepG2", 
                   forkPosition = 5,
                   flipMatrix = FALSE)
plotFPWM(fpwm,pdfName="MM1_HSA_HepG2_JUND_FPWM.pdf")
#
## Writing FPWM
#Save the FPWM to a local file can be used in matrix scanning or matrix clustering in transfact format or FPWMtransfact format. Transfact format will have a matrix for each interacting partner in the FPWM, while FPWMtransfact will output a single matrix.
write.FPWM(FPWM = fpwm, format = "transfac", fileName = "MM1_HSA_HepG2_JUND_FPWM.transfact" )
write.FPWM(FPWM = fpwm, format = "FPWMtransfac", fileName = "MM1_HSA_HepG2_JUND_FPWM.FPWMtransfac" )
#
#
## ForkedTF Case study 2: CEBPB in K562 cells
# Create a subdirectory for storing the output of the K562 CEBPB analysis.
sub_dir <- "K562_CEBPB"
dir.create(file.path(main_dir, sub_dir))
setwd(file.path(main_dir, sub_dir))
#
## MiniCofactor Report
# This function helps in the exploration of binding partners in a cell line. Input the mainTF, cell line of interest and the cobinding_threshold to generate a PDF report of the binding partners.
library(forkedTF)
# miniCofactorReport( TF = "CEBPB", cell = "K562" ) 
# Let's directly compute the -log10(Adjusted P-value) for the factors with the highest peak overlaps, using the enrichment test implemented in https://github.com/remap-cisreg/ReMapEnrich
miniCofactorReport(TF = "CEBPB",cell = "K562", filterBy="q.significance") 
# Execution takes several minutes...
#
## FPWM creation and plot
# Use the createFPWM function to extract the motif, from empirical datasets, that a TF uses when binding with a partner TF. plotFPWM helps in visualizing the FPWM.
fpwm <- createFPWM(mainTF ="CEBPB",
                   partners = c("ATF4","ATF7","ATF3","JUND","FOS","CEBPD"),
                   cell = "K562", 
                   forkPosition = 5,
                   flipMatrix = FALSE)
plotFPWM(fpwm,pdfName="MM1_HSA_K562_CEBPB_FPWM.pdf")
#
## Writing FPWM
#Save the FPWM to a local file can be used in matrix scanning or matrix clustering in transfact format or FPWMtransfact format. Transfact format will have a matrix for each interacting partner in the FPWM, while FPWMtransfact will output a single matrix.
write.FPWM(FPWM = fpwm, format = "transfac", fileName = "MM1_HSA_K562_CEBPB_FPWM.transfact" )
write.FPWM(FPWM = fpwm, format = "FPWMtransfac", fileName = "MM1_HSA_K562_CEBPB_FPWM.FPWMtransfac" )
#
## Cite
# A manuscript describing forkedTF has been submitted. If you are currently using forkedTF, please cite us as follows:
# Tirado-Magallanes R, Ghayour-Khiavi A, Santana-Garcia W, Dyer M, Lin QXX, Usefi H, Jha S, Thomas-Chollier M, Thieffry D, Benoukraf T. Representing Transcription Factor Dimer Binding Sites Using Forked-Position Weight Matrices and Forked-Sequence Logos. version: 1.2.0, 2023, [website: https://github.com/benoukraflab/forkedTF]
#
## License
# This project is licensed under GNU General Public License - see LICENSE.txt for details.
