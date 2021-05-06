print.results <- function(data, name){
  
  temp <- data 
  
  filename = paste0("results/",name,".txt")
  sink(file=filename)
  print(temp)
  sink()
  print("Done writing to file")
  
  
}

CleanData <- function(phy, data) {
  original_taxa <- row.names(phy)
  if_diffs <- name.check(phy, data)
  
  treedata(phy, data, sort = TRUE, warnings = TRUE)
}  

VisualizeData <- function(phy, data, file_1, file_2) {
  
  plot(phy, type = "fan", show.tip.label=TRUE, edge.width = 0.1)
  newdata <- data
  print(newdata)
  
  plot_tree(phy, file_1, TRUE)
  write.csv(newdata, file = file_2)
  
}

plot_tree <- function(file, tree, type, label, cex) {
  pdf(file=file)
  plot(tree, show.tip.label=label, cex=cex)
  dev.off()
}

convert.test.data <- function(data.table){
  
  temp <- data.table
  number.1 <- 0
  number.2 <- 0
  number.3 <- 0
  
  for(j in 1:dim(temp)[1]){
    
    for(i in 1:(dim(temp)[2]-1)){
      
      if(temp[[j,i]] == 0){
        
        if(temp[[j,i+1]] > 0){
          
          print("Condition 1")
          number.1 <- number.1 + 1
          temp[[j,i]] <- 1
          
        }
        
        else{
          
          print("Condition 3")
          number.3 <- number.3 +1
          temp[[j,i]] <- 3
        }
      }
      
      else{
        
        challenge <- (temp[[j,i+1]] - temp[[j,i]])/temp[[j,i]]
        
        if(challenge > 0.05){
          
          print("Condition 1")
          number.1 <- number.1 +1
          temp[[j,i]] <- 1
          
        }
        
        else if(challenge < -0.05){
          
          print("Condition 2")
          number.2 <- number.2 +1
          temp[[j,i]] <- 2
          
        }
        
        else{
          
          print("Condition 3")
          number.3 <- number.3 +1
          temp[[j,i]] <- 3
        }
        
      }
      
    }
    
  }
  
  sum.all <- number.1 + number.2 + number.3
  print(paste(number.1, number.2, number.3, sum.all))
  
  temp <- temp[,-dim(temp)[2]]
  colnames(temp) <- c("Value.1", "Value.2", "Value.3", "Value.4", "Value.5", 
                     "Value.6", "Value.7", "Value.8")
  return(temp)
  
}

assign.states <- function(data){
  
  temp <- data
  
  is.11 <- 0
  is.12 <- 0
  is.13 <- 0
  is.21 <- 0
  is.22 <- 0
  is.23 <- 0
  is.31 <- 0
  is.32 <- 0
  is.33 <- 0
  
  for(i in 1:dim(temp)[1]){
    
    state<-paste0(temp[i,1],temp[i,2])
    
    temp[[i,1]]<- c(state)
    
    if(temp[[i,1]]==11){
      is.11 <- 1
    }
    else if(temp[[i,1]]==12){
      is.12 <-1
    }
    else if(temp[[i,1]]==13){
      is.13 <-1
    }
    else if(temp[[i,1]]==21){
      is.21 <-1
    }
    else if(temp[[i,1]]==22){
      is.22 <-1
    }
    else if(temp[[i,1]]==23){
      is.23 <-1
    }
    else if(temp[[i,1]]==31){
      is.31 <-1
    }
    else if(temp[[i,1]]==32){
      is.32 <-1
    }
    else if(temp[[i,1]]==33){
      is.33 <-1
    }
    else{
      
      print(paste("Unknown value -",temp[i,1]))
      
    }
  }
  
  not.exist <- vector()
  
  if(is.11==0){
    not.exist <- append(not.exist, 11, length(not.exist))
  }
  if(is.12==0){
    not.exist <- append(not.exist, 12, length(not.exist))
  }
  if(is.13==0){
    not.exist <- append(not.exist, 13, length(not.exist))
  }
  if(is.21==0){
    not.exist <- append(not.exist, 21, length(not.exist))
  }
  if(is.22==0){
    not.exist <- append(not.exist, 22, length(not.exist))
  }
  if(is.23==0){
    not.exist <- append(not.exist, 23, length(not.exist))
  }
  if(is.31==0){
    not.exist <- append(not.exist, 31, length(not.exist))
  }
  if(is.32==0){
    not.exist <- append(not.exist, 32, length(not.exist))
  }
  if(is.33==0){
    not.exist <- append(not.exist, 33, length(not.exist))
  }
  
  temp2 <- as.matrix(temp[,1])
  row.names(temp2) <- row.names(temp)
  colnames(temp2) <- c("states")
  
  both <- list(state.data = temp2, missing = not.exist)
  return(both)
}

assign.control <- function(data){
  
  temp <- data
  temp2 <- as.matrix(temp[,2])
  row.names(temp2) <- temp[,1]
  colnames(temp2) <- c("states")
  
  return(temp2)
  print(temp2)
  
}

plot_sim_tree <- function(simm, label) {
  filename = paste0("results/SimmMap_",label,".pdf")
  pdf(file=filename)
  phytools::plotSimmap(simm[[1]], fsize = 0.5)
  dev.off()
}

print_plotAnc <- function(phyDat, recon, type,cex) {
  
    phangorn::plotAnc(phyDat, recon,1, cex.pie = par("cex.sub"), cex.sub = 0.4, cex=cex)
    
    filename = paste0("results/",type,".pdf")
    pdf(file=filename)
    phangorn::plotAnc(phyDat, recon,1, cex.pie = par("cex.sub"), cex.sub = 0.4, cex=cex)
    dev.off()

}

print_rates <- function(rates, variable) {
  
  filename = paste0("results/transition_rates_",variable)
  sink(file=filename)
  print(rates)
  sink(file=NULL)
  print("Done writing to file")
  
}

replace_NA <- function(blah) {
  
  model = blah
  temp_matrix <- matrix(data = NA, nrow = 2, ncol = 2)
  temp_matrix[[1,2]] <- model[[1,2]]
  print(temp_matrix[[1,2]])
  temp_matrix[[2,1]] <- model[[2,1]]
  print(temp_matrix[[2,1]])
  temp_matrix[[1,1]] <- 0
  temp_matrix[[2,2]] <- 0
  print(temp_matrix)
}

convert.test.2 <- function(test_matrix, ncol, name){
  
  temp <- test_matrix
  final <- matrix(NA, nrow = length(temp)+1, ncol = ncol+1)
  
  final[1,] <- name
  
  for(i in 1:length(temp)){
    
    current <-i
    new <- i+1
    
    final[[new,1]] <- paste0(names(temp)[current],"-2")
    final[[new,2]] <- temp[[current]][,1] + temp[[current]][,5] + temp[[current]][,6] + temp[[current]][,7]
    final[[new,3]] <- temp[[current]][,3]
    final[[new,4]] <- temp[[current]][,8]
    final[[new,5]] <- temp[[current]][,9]
    final[[new,6]] <- temp[[current]][,4]
    final[[new,7]] <- temp[[current]][,2]
  
  }
  
  final <- t(final)
  colnames(final) <- as.matrix(final[1,])
  final <- final[-1,]
  row.names(final) <- as.matrix(final[,1])
  final <- final[,-1]
  
  class(final) <- "numeric"
  
  return(final)
}


convert.control <- function(control, ncol, name){
  
  temp <- control
  final <- matrix(NA, nrow = length(temp)+1, ncol = ncol+1)
  
  final[1,] <- name
  
  for(i in 1:length(temp)){
    
    final[[i+1,1]] <- names(temp)[i]
    final[[i+1,2]] <- temp[[i]][,1] 
    final[[i+1,3]] <- temp[[i]][,2]
    final[[i+1,4]] <- temp[[i]][,3]
    final[[i+1,5]] <- temp[[i]][,4]
    final[[i+1,6]] <- temp[[i]][,5]
    final[[i+1,7]] <- temp[[i]][,6]
    
  }
  
  final <- t(final)
  colnames(final) <- as.matrix(final[1,])
  final <- final[-1,]
  row.names(final) <- as.matrix(final[,1])
  final <- final[,-1]
  
  class(final) <- "numeric"
  
  return(final)
}

convert_to_OTU_table <- function(control, testing){
  
  con <- control
  test <- testing
  
  con <- cbind(con, test)
  
  print(con)
  return(con)
}

create_sample_matrix <- function(OTU, name1, name2){
  
  temp <- OTU
  blah <- matrix(NA, nrow = ncol(OTU), ncol = 3)
  
  colnames(blah) <- c("method","points","bacteria")
  row.names(blah) <- t(colnames(temp))
  
  for(i in 1:(nrow(blah)/2)){
    
    blah[[i,1]] <- name1
    blah[[i,2]] <- 2
    blah[[i,3]] <- colnames(temp)[i]
  }
  
  for(i in ((nrow(blah)/2)+1):nrow(blah)){
    
    blah[[i,1]] <- name2
    blah[[i,2]] <- 2
    blah[[i,3]] <- colnames(temp)[i-(nrow(blah)/2)]
  }
  
  return(blah)
}

convert.pick.2 <- function(pick_matrix, ncol, name){
  
  temp <- pick_matrix
  final <- matrix(NA, nrow = length(temp)+1, ncol = ncol+1)
  
  final[1,] <- name
  
  for(i in 1:length(temp)){
    
    current <-i
    new <- i+1
    
    final[[new,1]] <- paste0(names(temp)[current],"-3")
    final[[new,2]] <- temp[[current]][,1] + temp[[current]][,2] + temp[[current]][,3] + temp[[current]][,8]
    final[[new,3]] <- temp[[current]][,5]
    final[[new,4]] <- temp[[current]][,6]
    final[[new,5]] <- temp[[current]][,9]
    final[[new,6]] <- temp[[current]][,4]
    final[[new,7]] <- temp[[current]][,7]
    
  }
  
  final <- t(final)
  colnames(final) <- as.matrix(final[1,])
  final <- final[-1,]
  row.names(final) <- as.matrix(final[,1])
  final <- final[,-1]
  
  class(final) <- "numeric"
  
  return(final)
}

create_sample_matrix_all <- function(OTU, name1, name2, name3){
  
  temp <- OTU
  blah <- matrix(NA, nrow = ncol(OTU), ncol = 3)
  
  colnames(blah) <- c("method","points","bacteria")
  row.names(blah) <- t(colnames(temp))
  
  stop1 <- nrow(blah)/3
  print(stop1)
  
  stop2 <- (nrow(blah)/3)*2
  print(stop2)
  
  for(i in 1:stop1){
    
    blah[[i,1]] <- name1
    blah[[i,2]] <- 2
    blah[[i,3]] <- colnames(temp)[i]
  }
  
  for(i in (stop1+1):stop2){
    
    blah[[i,1]] <- name2
    blah[[i,2]] <- 2
    blah[[i,3]] <- colnames(temp)[i-stop1]
  }
  
  for(i in (stop2+1):nrow(blah)){
    
    blah[[i,1]] <- name3
    blah[[i,2]] <- 2
    blah[[i,3]] <- colnames(temp)[i-stop2]
  }
  
  return(blah)
}

assign.states.3 <- function(data){
  
  temp <- data
  
  is.111 <- 0
  is.112 <- 0
  is.113 <- 0
  is.121 <- 0
  is.122 <- 0
  is.123 <- 0
  is.131 <- 0
  is.132 <- 0
  is.133 <- 0
  is.211 <- 0
  is.212 <- 0
  is.213 <- 0
  is.221 <- 0
  is.222 <- 0
  is.223 <- 0
  is.231 <- 0
  is.232 <- 0
  is.233 <- 0
  is.311 <- 0
  is.312 <- 0
  is.313 <- 0
  is.321 <- 0
  is.322 <- 0
  is.323 <- 0
  is.331 <- 0
  is.332 <- 0
  is.333 <- 0
  
  
  for(i in 1:dim(temp)[1]){
    
    state<-paste0(temp[i,1],temp[i,2],temp[i,3])
    
    temp[[i,1]]<- c(state)
    
    if(temp[[i,1]]==111){
      is.111 <- 1
    }
    else if(temp[[i,1]]==112){
      is.112 <-1
    }
    else if(temp[[i,1]]==113){
      is.113 <-1
    }
    else if(temp[[i,1]]==121){
      is.121 <-1
    }
    else if(temp[[i,1]]==122){
      is.122 <-1
    }
    else if(temp[[i,1]]==123){
      is.123 <-1
    }
    else if(temp[[i,1]]==131){
      is.131 <-1
    }
    else if(temp[[i,1]]==132){
      is.132 <-1
    }
    else if(temp[[i,1]]==133){
      is.133 <-1
    }
    else if(temp[[i,1]]==211){
      is.211 <-1
    }
    else if(temp[[i,1]]==212){
      is.212 <-1
    }
    else if(temp[[i,1]]==213){
      is.213 <-1
    }
    else if(temp[[i,1]]==221){
      is.221 <-1
    }
    else if(temp[[i,1]]==222){
      is.222 <-1
    }
    else if(temp[[i,1]]==223){
      is.223 <-1
    }
    else if(temp[[i,1]]==231){
      is.231 <-1
    }
    else if(temp[[i,1]]==232){
      is.232 <-1
    }
    else if(temp[[i,1]]==233){
      is.233 <-1
    }
    else if(temp[[i,1]]==311){
      is.311 <-1
    }
    else if(temp[[i,1]]==312){
      is.312 <-1
    }
    else if(temp[[i,1]]==313){
      is.313 <-1
    }
    else if(temp[[i,1]]==321){
      is.321 <-1
    }
    else if(temp[[i,1]]==322){
      is.322 <-1
    }
    else if(temp[[i,1]]==323){
      is.323 <-1
    }
    else if(temp[[i,1]]==331){
      is.331 <-1
    }
    else if(temp[[i,1]]==332){
      is.332 <-1
    }
    else if(temp[[i,1]]==333){
      is.333 <-1
    }
    else{
      
      print(paste("Unknown value -",temp[i,1]))
      
    }
  }
  
  not.exist <- vector()
  
  if(is.111==0){
    not.exist <- append(not.exist, 111, length(not.exist))
  }
  if(is.112==0){
    not.exist <- append(not.exist, 112, length(not.exist))
  }
  if(is.113==0){
    not.exist <- append(not.exist, 113, length(not.exist))
  }
  if(is.121==0){
    not.exist <- append(not.exist, 121, length(not.exist))
  }
  if(is.122==0){
    not.exist <- append(not.exist, 122, length(not.exist))
  }
  if(is.123==0){
    not.exist <- append(not.exist, 123, length(not.exist))
  }
  if(is.131==0){
    not.exist <- append(not.exist, 131, length(not.exist))
  }
  if(is.132==0){
    not.exist <- append(not.exist, 132, length(not.exist))
  }
  if(is.133==0){
    not.exist <- append(not.exist, 133, length(not.exist))
  }
  if(is.211==0){
    not.exist <- append(not.exist, 211, length(not.exist))
  }
  if(is.212==0){
    not.exist <- append(not.exist, 212, length(not.exist))
  }
  if(is.213==0){
    not.exist <- append(not.exist, 213, length(not.exist))
  }
  if(is.221==0){
    not.exist <- append(not.exist, 221, length(not.exist))
  }
  if(is.222==0){
    not.exist <- append(not.exist, 222, length(not.exist))
  }
  if(is.223==0){
    not.exist <- append(not.exist, 223, length(not.exist))
  }
  if(is.231==0){
    not.exist <- append(not.exist, 231, length(not.exist))
  }
  if(is.232==0){
    not.exist <- append(not.exist, 232, length(not.exist))
  }
  if(is.233==0){
    not.exist <- append(not.exist, 233, length(not.exist))
  }
  if(is.311==0){
    not.exist <- append(not.exist, 311, length(not.exist))
  }
  if(is.312==0){
    not.exist <- append(not.exist, 312, length(not.exist))
  }
  if(is.313==0){
    not.exist <- append(not.exist, 313, length(not.exist))
  }
  if(is.321==0){
    not.exist <- append(not.exist, 321, length(not.exist))
  }
  if(is.322==0){
    not.exist <- append(not.exist, 322, length(not.exist))
  }
  if(is.323==0){
    not.exist <- append(not.exist, 323, length(not.exist))
  }
  if(is.331==0){
    not.exist <- append(not.exist, 331, length(not.exist))
  }
  if(is.332==0){
    not.exist <- append(not.exist, 332, length(not.exist))
  }
  if(is.333==0){
    not.exist <- append(not.exist, 333, length(not.exist))
  }
  
  temp2 <- as.matrix(temp[,1])
  row.names(temp2) <- row.names(temp)
  colnames(temp2) <- c("states")
  
  both <- list(state.data = temp2, missing = not.exist)
  return(both)
}

convert.test.3 <- function(test_matrix, ncol, name){
  
  temp <- test_matrix
  final <- matrix(NA, nrow = length(temp)+1, ncol = ncol+1)
  
  final[1,] <- name
  
  for(i in 1:length(temp)){
    
    current <-i
    new <- i+1
    
    final[[new,1]] <- paste0(names(temp)[current],"-triple")
    final[[new,2]] <- (temp[[current]][,1] + temp[[current]][,2] + temp[[current]][,3] 
                       + temp[[current]][,7] + temp[[current]][,12] + temp[[current]][,13])
    final[[new,3]] <- temp[[current]][,5] + temp[[current]][,8]
    final[[new,4]] <- temp[[current]][,14] + temp[[current]][,15]
    final[[new,5]] <- temp[[current]][,16] + temp[[current]][,17]
    final[[new,6]] <- temp[[current]][,4] + temp[[current]][,6] + temp[[current]][,9]
    final[[new,7]] <- temp[[current]][,10] + temp[[current]][,11]
    
  }
  
  final <- t(final)
  colnames(final) <- as.matrix(final[1,])
  final <- final[-1,]
  row.names(final) <- as.matrix(final[,1])
  final <- final[,-1]
  
  class(final) <- "numeric"
  
  return(final)
}