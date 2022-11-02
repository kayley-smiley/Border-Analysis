#the F function for the Elliptic Scan Method


ba_ellipse <- function(scan_output, nsim = 1000, coords = coords, cases = cases, pop = pop,
                      ex = sum(cases)/sum(pop) * pop, longlat = FALSE){

  coords <- as.matrix(coords)
  N <- nrow(coords)
  d <- sp::spDists(coords, longlat = longlat)
  ty <- sum(cases)
  
  nangle <- scan_output$rel_param$nangles
  shape <- scan_output$rel_param$shapes
  a <- scan_output$rel_param$a_penalty
  
  enn <- smerc:::elliptic.nn(coords, pop, ubpop = scan_output$rel_param$ubpop, shape, nangle)
  pen <- elliptic.penalty(a = a, enn$shape_nn) #change to smerc:::elliptic.penalty
 
  #define ein
  ein <- nn.cumsum(enn$nn, ex, simplify = FALSE)
  
  #compute eout
  eout <- lapply(ein, function(i){
    sum(ex) - unlist(i)
  })
  
  sim_res <- pbapply::pblapply(seq_len(nsim), function(i){
    
    #generate a multinomial data set denoted by ysim 
    ysim = stats::rmultinom(1, size = ty, prob = cases)
    
    #compute yin (keep in nn format)
    yin = nn.cumsum(enn$nn, ysim, simplify = FALSE)
    
    #compute t_obs (keep in nn format)
    tobs <- lapply(seq(1, length(enn$nn)), function(i){
      stat.poisson.adj(yin[[i]], ty, log(ein[[i]]), log(eout[[i]]), a = a,
                       pen = pen[i], min.cases = 2, return.max = FALSE)
    })
    #change to smerc:::stat.poisson.adj
    
    clusters <- noc_nn(enn$nn, tobs, nnoc = length(scan_output$clusters)) #change to smerc:::noc_nn
     
    return(unlist(clusters$clusts))
    
  })
  
  #find f_function values
  counts <- c(unlist(sim_res), seq(1,length(pop)))
  table_counts <- as.data.frame((table(counts) -1)/nsim)
  f_values <- table_counts$Freq
  
  return(f_values)
}

