rm(list = ls())
root = "S:/SHARE/drf/Automated_Nowcast/rujun/Code/new_code/Replicate/"
est_dir = "S:/SHARE/drf/Automated_Nowcast/rujun/Code/new_code/Replicate/estimation"
setwd(root)

library("Haver")
library("solidearthtide")
library("zoo")
library("R.matlab")
library("solidearthtide")
source("ts_trans.R")
library(lubridate)
library("signal")
library(MASS)
library(magic)

source("News_DFM.R")
source("para_const.R")
source("runKF.R")
source("DFM.R")
source("utils.R")

path_old = paste(root, "Old_Data", sep = "")
setwd(path_old)

load("Data_2016_11_01.RData")

data = Data_2016_11_01$data

setwd(root)

#data = data[, 1:34]
blocks = readMat("blocks.mat")
blocks = blocks$blocks
#blocks = blocks[1:34,]


Par = list()
# specify quarter data
Par$nQ <- 3
# specify local blocks
Par$blocks <- blocks
# specify number of factors to extract for each block
Par$r <- c(1, 1, 1, 1)
# specify number of lags
Par$p <- 1

res = DFM(data, Par)

