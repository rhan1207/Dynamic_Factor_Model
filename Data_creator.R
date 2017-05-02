############################################################################### 
##### Fetch data from Haver DLX and save new data and info in .RData file #####
############################################################################### 

startdate = "2000-01-01"
start_data = "2000-10-01"


haver.path()


lgd = vector()
Data = c()
legend = list()

for (i in 1:length(Vars)){
  x = dlxGetInfo(Vars[i], Extensions[i])
  lgd[i] = x$frequency
  legend_temp = list(x$frequency, x$dateTimeMod, units[[i]], tran[[i]], x$varname, x$endDate)
  names(legend_temp) = c("frequency", "dtmod", "units", "tran", "mnemonic", "period")
  legend = append(legend, list(legend_temp))
  d = dlxGetData(Vars[i], Extensions[i], startdate, today)
  temp = as.vector(coredata(d))

  if (x$frequency == 30){
    
    temp_idx = 3*c(1:length(temp))
    D = matrix(NA, length(temp) * 3, 1)
    D[temp_idx] = temp
    Data = cbind(Data, D[1:dim(Data)[1]])
  }
  else{
    Data = cbind(Data, temp)
  }
}

#transformed data
Data_new = matrix(NA, dim(Data)[1]-1, dim(Data)[2])
for (i in 1:dim(Data)[2]){
  if (units[i] == "lvl"){
    Data_new[,i] = lvl(Data[,i])
  }
  else if(units[i] == "pc"){
    Data_new[,i] = pc(Data[,i])
  } 
  else if(units[i] == "pca"){
    Data_new[,i] = pca(Data[,i])
  } 
  else if(units[i] == "pcq"){
    Data_new[,i] = pcq(Data[,i])
  } 
  else if(units[i] == "pcqa"){
    Data_new[,i] = pcqa(Data[,i])
  } 
  else if(units[i] == "lcpp"){
    Data_new[,i] = lcpp(Data[,i])
  } 
  
}


#identifying quarterly variables
indx=which(lgd == 30)

j = 1
Header = header
for (i in append(setdiff(1:length(lgd),indx), indx)){
  Header[j] = header[i]                                   # ordered Header, with quarterly vars at the end
  j = j + 1
}

Data_new = Data_new[, append(setdiff(1:length(lgd),indx), indx)]

start_idx = elapsed_months(start_data, startdate)

Data_new = Data_new[start_idx:dim(Data_new)[1], ]
# delete the last row if all are nan
if (sum(is.na(Data_new[dim(Data_new)[1],])) == Nseries){
  Data_new = Data_new[1:(dim(Data_new)[1]-1),]
}

### Saving data and legend
filename = paste("Data", substr(today, 1,4), substr(today, 6,7), substr(today, 9,10), sep="_")

data = list()
data$data = Data_new
data$legend = legend
assign(filename, data)

setwd(path_old)
save(list = filename, file = paste(filename, ".RData", sep=""))



