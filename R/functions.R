print.results <- function(data, name){
  
  temp <- data 
  
  filename = paste0("results/",name)
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
  filename = paste0("results/SimmMap_",label)
  pdf(file=filename)
  phytools::plotSimmap(simm[[1]], fsize = 0.5)
  dev.off()
}

print_plotAnc <- function(phyDat, recon, type,cex) {
  
    phangorn::plotAnc(phyDat, recon,1, cex.pie = par("cex.sub"), cex.sub = 0.4, cex=cex)
    
    filename = paste0("results/",type)
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