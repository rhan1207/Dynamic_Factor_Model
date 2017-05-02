# funtion defined in this file is only useful for vectors


lvl <- function(x){
  y = x[2:length(x)]
  return(y)
}

pc <- function(x){
  
  y = matrix(0, length(x)-1, 1)
  for (i in 2:length(x)){
    y[i-1] = 100*(x[i] - x[i-1])/x[i-1]
  }
  return(y)
}

pca <- function(x){
  
  y = matrix(0, length(x)-1, 1)
  for (i in 4:length(x)){
    y[i-1] = 400*(x[i] - x[i-3])/x[i-3]
  }
  return(y)
}

pcq <- function(x){
  
  y = matrix(0, length(x)-1, 1)
  for (i in 4:length(x)){
    y[i-1] = 100*(x[i] - x[i-3])/x[i-3]
  }
  return(y)
}

pcqa <- function(x){
  
  y = matrix(0, length(x)-1, 1)
  for (i in 4:length(x)){
    y[i-1] = 400*(x[i] - x[i-3])/x[i-3]
  }
  return(y)
}

lcpp <- function(x){
  
  y = matrix(0, length(x)-1, 1)
  for (i in 2:length(x)){
    y[i-1] = x[i] - x[i-1]
  }
  return(y)
}


elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}