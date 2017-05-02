rm(list = ls())
root = "S:/SHARE/drf/Automated_Nowcast/rujun/Code/new_code/Replicate/"
est_dir = "S:/SHARE/drf/Automated_Nowcast/rujun/Code/new_code/Replicate/estimation"

library(ggplot2)

path_old = paste(root, "Old_Data", sep = "")
setwd(path_old)

load("Data_2016_11_01.RData")

data = Data_2016_11_01$data

data = data[, 1:34]

setwd(root)
load("1.RData")
load("2.RData")
load("3.RData")

temp_df = data.frame(data[,1], aaa[,1], bbb[,1], ccc[,1])

plot(1:193, aaa[,2], type = "p")
lines(1:193, bbb[,2], type = "p", col = "red")