#the F function for the Circular Scan Method

ba_circle <- function(scan_output, nsim = 1000, coords = coords, cases = cases, pop = pop,
                      ex = sum(cases)/sum(pop) * pop, longlat = FALSE){
  if(scan_output$clusters[[1]]$pvalue > 0.05){
    f_values <- rep(0, length(cases))
  } else{
  coords <- as.matrix(coords)
  N <- nrow(coords)
  d <- sp::spDists(coords, longlat = longlat)
  nn <- smerc::scan.nn(d, pop, scan_output$rel_param$ubpop)
  ty <- sum(cases)
  
  #define ein
  ein <- smerc::nn.cumsum(nn, ex, simplify = FALSE)
  
  #compute eout
  eout <- lapply(ein, function(i){
    sum(ex) - unlist(i)
  })
  
  #simulations & results
  sim_res <- pbapply::pblapply(seq_len(nsim), function(i){
    
    #generate a multinomial data set denoted by ysim 
    ysim = stats::rmultinom(1, size = ty, prob = cases)
    
    #compute yin (keep in nn format)
    yin = smerc::nn.cumsum(nn, ysim, simplify = FALSE)
    
    #compute t_obs (keep in nn format)
    tobs <- lapply(seq(1, length(ex)), function(i){
      stat.poisson(yin[[i]], ty - yin[[i]], ein[[i]], eout[[i]])
    })
    
    clusters <- noc_nn(nn, tobs, nnoc = length(scan_output$clusters)) #noc_nn in smerc
    
    #return clusters
    return(unlist(clusters$clusts))
  })
  
  #calculate f_function values
  counts <- c(unlist(sim_res), seq(1,length(nn)))
  table_counts <- as.data.frame((table(counts) -1)/nsim)
  f_values <- table_counts$Freq
  }
  
  return(f_values)
}
