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

source("News_DFM.R")
source("para_const.R")
source("runKF.R")

path_old = paste(root, "Old_Data", sep = "/")

### Organizing data sets by date ###

#checking date of oldest data set we used for nowcast
old_files = list.files(path = path_old)
tmp = old_files[length(old_files)]

yesterday = as.Date(paste(substr(tmp, 6, 9), substr(tmp, 11, 12), substr(tmp, 14, 15), sep="/"))
old_date = Datenum(yesterday)

# today's date
today = Sys.Date() # use this date string in Data_creator
new_date =  Datenum(today)

#maimum forward forecast horizon
extra_dates = 10


### Importing Legend, header and preparing for reordering ###

# Loading current specification
hd = readMat("header.mat")$header
header = c()
for (h in 1:length(hd)){
  header = append(header, hd[[h]][[1]])
}


Spec = read.csv("spec2.csv", header=FALSE, sep=",", strip.white=TRUE)

Vars = lapply(as.vector(Spec[,1]), function(x) substr(x, 2, nchar(x)-1))
Extensions = lapply(as.vector(Spec[,2]), function(x) substr(x, 2, nchar(x)-1))
tran = lapply(as.vector(Spec[,4]), function(x) substr(x, 2, nchar(x)-1))
units =  lapply(as.vector(Spec[,5]), function(x) substr(x, 2, nchar(x)-1))


Nseries = length(header)

if (old_date != new_date){
  print("Updating data...")
  source("Data_creator.R")
  print("update finished.")
  
  Legend = legend
  Data_New = data$data
  load(old_files[length(old_files)])
  eval(parse(text = paste("Data_Old =", substr(old_files[length(old_files)], 1, 15))))
  Data_Old = Data_Old$data
  setwd(root)
}else{
  print("No update needed.")
  quit()
}


#checking if new data is actually a change w/ respect to old data
if (identical(Data_New, Data_Old)){
  print("No data series were updated")
  setwd(path_old)
  file.remove(paste(filename, ".RData", sep=""))
  quit()
}


# adding lines to match sizes of Data_Old and Data_New
r_old = dim(Data_Old)[1]
c_old = dim(Data_Old)[2]

r_new = dim(Data_New)[1]
c_new = dim(Data_New)[2]


while ((r_new - r_old) > 0){
  Data_Old = rbind(Data_Old, rep(NA, Nseries))
  r_old = r_old + 1
}

### Reordering: placing qurterly series at the end of the data set ###

use_tobe = c(1:Nseries)
use_tobe = t(rbind(use_tobe, append(setdiff(c(1:length(lgd)), indx), indx)))

#position of GDP among vars
v_news = which(Header == "Real gross domestic product (saar, bil.chn.2009$)")


#Estimation
if (as.yearqtr(yesterday) != as.yearqtr(today)){
  estim_choice = 1
  setwd(est_dir)                  
  Res_Old = load("Res.RData")
  save(Res_Old, "Res_Old.RData")
  print("Starting estimation...")
  blocks = readMat("blocks.mat")
  blocks = blocks$blocks
  

  Par = list()
  Par$nQ <- 3
  Par$blocks <- blocks
  Par$r <- c(1, 1, 1, 1)
  Par$p <- 1
  browser()
  Data_est = Data_New[(dim(data_new)[1] - 12*15 + 1):dim(data_new)[1],]
  Res = DFM(Data_est, Par)
  
  save(Res, file = "Res.RData")
  
  est_archive = paste(root, "estimation", "Archive", sep = "/")
  setwd(est_archive)
  save(Res, file = paste(substr(as.yearqtr(today), 1,4), substr(as.yearqtr(today), 6,7), sep = "/"))

  setwd(root)
  print("Estimation ending.")
}else{
  estim_choice = 0
  setwd(est_dir)
  load("Res.RData")
}

setwd(root)
### Adding nans to data in order to forecast
i = 1 
while (i <= extra_dates){
  Data_New = rbind(Data_New, matrix(NA, 1, Nseries))
  Data_Old = rbind(Data_Old, matrix(NA, 1, Nseries))
  i = i + 1
}

### Data Revisions and News for prevoius quarter, current quarter and next quarter

curr_year = strtoi(substring(as.yearqtr(today), 1,4))
curr_quart = strtoi(substring(as.yearqtr(today), 7, 7))

DatesV = seq(as.Date(startdate), by = "month", length.out = (dim(Data_New)[1] + 20))
DatesV = DatesV[10:length(DatesV)]



for (i in (-1:1)){       
  if (i==-1 & curr_quart==1){
    quart = 4
    yr = curr_year - 1
  }else{
    if ((curr_quart + i) %% 4 == 0){
      quart = 4
      yr = curr_year
    }else{
      quart = (curr_quart + i) %% 4
      yr = curr_year + floor((curr_quart+i)/4)
    }
  }
  file_name = paste(toString(yr), "Q", toString(quart), ".csv", sep = "")
  
  date_fcst = Datenum(as.Date(paste((toString(yr)), toString(quart * 3), "01", sep="/")))
  
  t_fcst = match(date_fcst, Datenum(DatesV))
  
  
  # check if GDP release has happened
  #A1 = Data_New[which(!is.na(Data_New[, v_news])), v_news]
  #A2 = Data_Old[which(!is.na(Data_New[, v_news])), v_news]
  
  
  load("GDP_cal.RData")
  whr_gdp = which(GDP_releases$quarter == substring(file_name, 1, 6))
  whr_gdp_tmp = which(GDP_releases$datenum[[1]][whr_gdp] <= new_date)
  
  GDP_today = 0
  
  if (length(whr_gdp_tmp) == 0){GDP_rel = 0}
  else{
    if (length(whr_gdp_tmp) == 1){GDP_rel = 1}
    else{ 
      if (length(whr_gdp_tmp) == 2){GDP_rel = 2}
      else{ 
        if (length(whr_gdp_tmp) == 3){GDP_rel = 3}
      }    
    }
  }
  
  ###
  if ((length(whr_gdp) > 0) & (GDP_rel > 0)){
    if (new_date == GDP_releases$datenum[[1]][whr_gdp[GDP_rel]]){GDP_today = 1}
  }
  
  
  if (GDP_rel == 0){
    
    # data revision
    Data_Old_rev = Data_New
    Data_Old_rev[is.na(Data_Old)] = NA
    
    res_rev = News_DFM(Data_Old_rev, Data_New, Res, t_fcst, v_news)
    y_old_rev = res_rev$y_old
    y_new = res_rev$y_new
    fore_rev = res_rev$fore
    actual_rev = res_rev$actual
    filt_rev = res_rev$filt
    
    #news
    res = News_DFM(Data_Old, Data_New, Res, t_fcst, v_news)
    y_old = res$y_old 
    
    data_rev = y_old_rev - y_old  
    if (estim_choice == 1){
      res_est = News_DFM(Data_Old, Data_New, Res_Old, t_fcst, v_news)
      y_old_est = res_est$y_old
      data_est = y_old - y_old_est
      y_old = y_old_est
    }
  }else{
    fore_rev = NA
    if ((GDP_today == 1) & (GDP_rel == 1)){
      res_GDP = News_DFM(Data_Old, Data_New, Res, t_fcst, v_news)
      y_new = res_GDP$y_new
      y_old = res_GDP$y_old
      data_rev = y_new - y_old 
    }
  }
  
  if ((GDP_rel == 0) | ((GDP_rel == 1) & (GDP_today == 1))){
  
  ### Display
  extra_row = 2 + estim_choice
  news_out = c()
  if (! is.na(fore_rev)){
    for (i in 1:Nseries){
      if (fore_rev[i,1] != 0) {
        news_out = cbind(news_out, i)
      }
    }
  }
  structure = data.frame(matrix(nrow = (length(news_out) + extra_row), ncol = 10))
  
  #Special structure for GDP release
  if (GDP_today == 1){
    if (GDP_rel == 1){structure = data.frame(matrix(nrow = 3, ncol = 10))}
    else{structure = data.frame(matrix(nrow = 2, ncol = 10))}
  }
  
  # output header
  structure[1,1] = paste(substring(Datevec(old_date), 6,7), substring(Datevec(old_date), 9, 10), substring(Datevec(old_date), 1, 4), sep = "/")
  structure[1,2] = ""
  structure[1,3] = ""
  structure[1,4] = ""
  structure[1,5] = ""
  structure[1,6] = ""
  structure[1,7] = ""
  structure[1,8] = ""
  structure[1,9] = ""
  structure[1,10] = y_old
  
  j = 0 # for GDP_today = 1 and GDP_rel > 1
  browser()
  # output content
  if ((GDP_rel == 0) | (GDP_today == 1)) {
  
    if (GDP_rel == 0){  
      for (j in 1:(length(news_out) + extra_row-1)){
        # data revision row
        if (j == (length(news_out) + 1)){
          structure[j+1,1] = ""
          structure[j+1,2] = ""
          structure[j+1,3] = "Data revisions"
          structure[j+1,4] = ""
          structure[j+1,5] = ""
          structure[j+1,6] = ""
          structure[j+1,7] = ""
          structure[j+1,8] = ""
          structure[j+1,9] = data_rev
          structure[j+1,10] = ""    

        } else {
        # parameter revision row
          if (j == (length(news_out) + 2)){
            structure[j+1,1] = ""
            structure[j+1,2] = ""
            structure[j+1,3] = "Parameter revisions"
            structure[j+1,4] = ""
            structure[j+1,5] = ""
            structure[j+1,6] = ""
            structure[j+1,7] = ""
            structure[j+1,8] = ""
            structure[j+1,9] = data_est
            structure[j+1,10] = "" 

          } else {
            k = news_out[j]
            structure[j+1,1] = ""
            structure[j+1,2] = paste(substr(Legend[[k]]$dtmod, 1,4), substr(Legend[[k]]$dtmod, 6,7), substr(Legend[[k]]$dtmod, 9,10), sep="/")
            structure[j+1,3] = Header[[k]]
            structure[j+1,4] = format(Legend[[k]]$period, "%b")
            structure[j+1,5] = Legend[[k]]$tran
            structure[j+1,6] = fore_rev[k,1]
            structure[j+1,7] = actual_rev[k,1]
            structure[j+1,8] = filt_rev[k,1,1]
            structure[j+1,9] = filt_rev[k,1,1] * (actual_rev[k,1]-fore_rev[k,1])
            structure[j+1,10] = ""
          }
        }
      }
    }else{
      # for first release, only GDP series matter;
      # else, the structure content is empty
      if (GDP_rel == 1){        
        k = v_news
        structure[2,1] = ""
        structure[2,2] = paste(substr(Legend[[k]]$dtmod, 1,4), substr(Legend[[k]]$dtmod, 6,7), substr(Legend[[k]]$dtmod, 9,10), sep="/")
        structure[2,3] = Header[[k]]
        structure[2,4] = format(Legend[[k]]$period, "%b")
        structure[2,5] = Legend[[k]]$tran
        structure[2,6] = y_old
        structure[2,7] = y_new
        structure[2,8] = 1
        structure[2,9] = data_rev
        structure[2,10] = ""  
        j = 1
      }
    }
  }
  # output ending
  structure[j+2,1] = paste(substring(Datevec(new_date), 6,7), substring(Datevec(new_date), 9, 10), substring(Datevec(new_date), 1, 4), sep = "/")
  structure[j+2,2] = ""
  structure[j+2,3] = ""
  structure[j+2,4] = ""
  structure[j+2,5] = ""
  structure[j+2,6] = ""
  structure[j+2,7] = ""
  structure[j+2,8] = ""
  structure[j+2,9] = ""
  structure[j+2,10] = y_new
  
  write.table(structure, file_name, col.names = FALSE, row.names=FALSE, append = T, sep = ",") 
  
  }
}
























