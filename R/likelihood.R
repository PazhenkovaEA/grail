# Likelihood for allele ratio. 
# Theoretical ratio is centered on m and with std (chosen by user). 
# We use a truncated normal distribution on [0,1. Expectation from LL for heterozygotes is substracted.

llratio_esp <- function(x, Ntot = 500, m = 0.5, sratio = 0.25){
  -(x-m)^2/2/(sratio^2+m/Ntot)+1/2-
    (m^2*exp(-m^2/2/(sratio^2+m/Ntot)) + (1-m)^2*exp(-(1-m)^2/2/(sratio^2+m/Ntot)))/
    sqrt((sratio^2+m/Ntot))/2/sqrt(2*pi)/     (pnorm(1, m, sqrt((sratio^2+m/Ntot)))-                                
                                                 pnorm(0, m, sqrt((sratio^2+m/Ntot))))
}

# The rate of reads of unexpected variants follows an Exponential distribution.
# The number of contamination reads depends on the genotype. 
# It takes into account any sequence that is not an estimated allele,
# a mutation at distance 1 from an allele or a slippage error at distance 2 or less from an allele.

llcontamination <- function(x, lambda = 1/0.05) log(lambda)-lambda*x


# Abundance of mutated reads ( for all sequences together.
# This quantity depends on a mutation rate via a Poisson distribution.

llabmut <- function(x, N1 = 250, mu = 0.1){log(dpois(x, mu*N1))}

# Number of different mutated sequences knowing the abundance of mutated reads

llvariantobs <- function(x, Nmut = 25){log(dpois(x, Nmut))}

# The number of reads obtained by slippage from A and B follows also a Poisson distribution.
# Sequence slippage (stutter) distance dA from A and dB from B.
# The calculation of this LL is done up to 2 of distance before A and after B (with A before B).
# Reads from A created by slipping from B (and vice versa) are ignored.

llslippage <- function(x, dA,
                       N1 = 250, N2 = 100, d12 = 3,
                       nu_minus = 0.05, nu_plus = 0.1){
  #dA : distance from A
  dB <- dA - d12
  if (dA < 0) nuA <- nu_minus
  else if (dA > 0) nuA <- nu_plus
  else {return(0)}
  if (dB < 0) nuB <- nu_minus
  else if (dB > 0) nuB <- nu_plus
  else {return(0)}
  log(dpois(x, nuA^abs(dA)*N1+nuB^abs(dB)*N2))
}


# calculate distance between alleles

dist12 <- function(data, genotype){
  if (nchar(data$sequences[genotype[1]]) < nchar(data$sequences[genotype[2]])){
    d12 <- 0
    graph <- genotype[1]
    while (!(genotype[2] %in% graph)){
      graph2 <- unique(c(graph, data$slips[data$slips$col %in% graph,]$row))
      d12 <- d12 + 1
      if (length(graph2) == length(graph)){
        print("link not found")
        return(NA)
      }
      graph <- graph2
    }
  }
  else {
    d12 <- 0
    graph <- genotype[1]
    while (!(genotype[2] %in% graph)){
      graph2 <- unique(c(graph, data$slips[data$slips$row %in% graph,]$col))
      d12 <- d12 - 1
      if (length(graph2) == length(graph)){
        print("link not found")
        return(NA)
      }
      graph <- graph2
    }
  }
  d12
}


#Looks for slipages

graph_slippage <- function(data, genotype, dslipmax, repli){
  d12 <- dist12(data, genotype)
  slippage1 <- data.frame(seq = genotype[1], dist1 = 0)
  for (d in 1:dslipmax){ #backward
    auxseq <- unique(data$slips[data$slips$row %in% slippage1[slippage1$dist1 == -(d-1),]$seq,]$col)
    if (length(auxseq) > 0) {
      slippage1 <- rbind(slippage1, data.frame(seq = auxseq, dist1 = -d))
    }
  }
  for (d in 1:dslipmax){ #forward
    auxseq <- unique(data$slips[data$slips$col %in% slippage1[slippage1$dist1 == (d-1),]$seq,]$row)
    if (length(auxseq) > 0) {
      slippage1 <- rbind(slippage1, data.frame(seq = auxseq, dist1 = d))
    }
  }
  if (genotype[1] != genotype[2]){
    slippage2 <- data.frame(seq = genotype[2], dist = 0)
    for (d in 1:dslipmax){ #backward
      auxseq <- unique(data$slips[data$slips$row %in% slippage2[slippage2$dist == -(d-1),]$seq,]$col)
      if (length(auxseq) > 0) {
        slippage2 <- rbind(slippage2, data.frame(seq = auxseq, dist = -d))
      }
    }
    for (d in 1:dslipmax){ #forward
      auxseq <- unique(data$slips[data$slips$col %in% slippage2[slippage2$dist == (d-1),]$seq,]$row)
      if (length(auxseq) > 0) {
        slippage2 <- rbind(slippage2, data.frame(seq = auxseq, dist = d))
      }
    }
    slippage2$dist1 <- slippage2$dist + d12
    slippage2$dist <- NULL
    
    slippage <- rbind(slippage1, slippage2)
  }
  else {
    slippage <- slippage1
  }
  
  slippage <- slippage[!duplicated(slippage),]
  slippage <- slippage[order(slippage$dist1),]
  rownames(slippage)<- NULL
  slippage$counts <- data$counts[slippage$seq,repli]
  
  slippage
}


#' Estimate genotype likelihood for each replicate separately.
#'
#' This function calculates genotype likelihood for all the replicates of the sample. 
#' @param data TODO: Describe the format
#' @param dslipmax Maximum number of motifs between allele and slipage.
#' @param repli Replicates to process.
#' @param filter Number of reads threshold.
#' @param param list. m = expected allele ratio for heterozygotes, sratio = standard deviation of the ratio for heterozygotes, contamination rate = expected proportion of contamination,
#' nu_plus = rate slippage +, nu_minus = rate: slippage -.
#' @keywords microsatellites genotype likelihood
#' @export
#' @examples
#' ll_sample()

lltotal <- function(data = testData,
                    genotype = c(1,1),
                    repli = 1,
                    dslipmax = 1,
                    filter = 0,
                    param = list(m = m,
                                 sratio = sratio,
                                 lambda = lambda,
                                 mu = mu,
                                 nu_plus = nu_plus,
                                 nu_minus = nu_minus)){
  
  Nreads <- sum(data$counts[,repli])
  if (Nreads <= filter){
    return(NULL)
  }
  heterozygote <- genotype[1] != genotype[2]
  
  if (heterozygote){
    d12 <-dist12(data, genotype)
    if (d12 < 0){
      d12 <- -d12
      genotype <- c(genotype[2], genotype[1])
    }
    
    N1 <- data$counts[genotype[1],repli]
    Var1 <- data$SNPs[data$SNPs$col == genotype[1],]$row
    Var1 <- Var1[Var1 != genotype[1] & data$counts[Var1,repli] != 0]
    Nmut1 <- sum(data$counts[Var1,repli])
    V1 <- length(Var1)
    
    N2 <- data$counts[genotype[2],repli]
    Var2 <- data$SNPs[data$SNPs$col == genotype[2],]$row
    Var2 <- Var2[Var2 != genotype[2] & data$counts[Var2,repli] != 0]
    Nmut2 <- sum(data$counts[Var2,repli])
    V2 <- length(Var2)
    
    d12 <-dist12(data, genotype)
  }
  else {
    N1 <- data$counts[genotype[1],repli]
    Var1 <- data$SNPs[data$SNPs$col == genotype[1],]$row
    Var1 <- Var1[Var1 != genotype[1] & data$counts[Var1,repli] != 0]
    Nmut1 <- sum(data$counts[Var1,repli])
    V1 <- length(Var1)
    
    d12 <- 0
  }
  
  slippage <- graph_slippage(data, genotype, dslipmax, repli)
  
  graph <- unique(c(genotype[1], Var1, slippage$seq))
  if (heterozygote){
    graph <- unique(c(graph, unique(c(genotype[2], Var2))))
  }
  
  ll <- NULL
  
  #ratio A/(A+B)
  if (heterozygote){
    aux <- llratio_esp(N1/(N1+N2),Ntot = N1+N2,m = param$m,sratio = param$sratio)
    if (is.nan(aux)) aux <- -Inf
    ll <- rbind(ll, data.frame(component = "ratio_1/(1+2)-expected",
                               ll = aux))
    
  }
  #contamination
  Nc <- Nreads - sum(data$counts[graph,repli])
  ll <- rbind(ll, data.frame(component = "contamination",
                             ll = llcontamination(Nc/Nreads, lambda = param$lambda)))
  #PCR errors
  if (heterozygote){
    ll <- rbind(ll, data.frame(component = "reads_mutations1",
                               ll = llabmut(Nmut1, N1 = N1, mu = param$mu)))
    ll <- rbind(ll, data.frame(component = "variants_mutations1",
                               ll = llvariantobs(V1, Nmut = Nmut1)))
    
    ll <- rbind(ll, data.frame(component = "reads_mutations2",
                               ll = llabmut(Nmut2, N1 = N2, mu = param$mu)))
    ll <- rbind(ll, data.frame(component = "variants_mutations2",
                               ll = llvariantobs(V2, Nmut = Nmut2)))
  }
  else {
    ll <- rbind(ll, data.frame(component = "2*reads_mutations1",
                               ll = 2*llabmut(Nmut1, N1 = N1, mu = param$mu)))
    ll <- rbind(ll, data.frame(component = "2*variants_mutations1",
                               ll = 2*llvariantobs(V1, Nmut = Nmut1)))
    
  }
  #slippage
  list_splippage <- unique(slippage$dist1)
  n_non0 <- length(list_splippage)-1
  if (heterozygote) n_non0 <- n_non0 - 1
  for (d1 in list_splippage){
    if (heterozygote){
      ll <- rbind(ll, data.frame(component = paste0("slippage_from1_", d1),
                                 ll = llslippage(sum(slippage[slippage$dist1 == d1,]$counts),
                                                 d1, N1 = N1, N2 = N2, d12 = d12)/n_non0))
    }
    else {
      ll <- rbind(ll, data.frame(component = paste0("slippage_from1_", d1),
                                 ll = llslippage(sum(slippage[slippage$dist1 == d1,]$counts),
                                                 d1, N1 = N1, N2 = 0, d12 = d12)/n_non0))
    }
  }
  ll <- rbind(data.frame(component = "total", ll = sum(ll$ll)), ll)
  ll
}

#' Estimate genotype likelihood for sample.
#'
#' This function calculates genotype likelihood for all the replicates of the sample. 
#' @param data TODO: Describe the format
#' @param dslipmax Maximum number of motifs between allele and slippage.
#' @param repli Replicates to process.
#' @param filter Number of reads threshold.
#' @param param list. m = expected allele ratio for heterozygotes, sratio = standard deviation of the ratio for heterozygotes, contamination rate = expected proportion of contamination,
#' nu_plus = rate slippage +, nu_minus = rate: slippage -.
#' @keywords microsatellites genotype likelihood
#' @export
#' @examples
#' ll_sample()


ll_sample <- function(data = testData,
                      genotype = c(1,1),
                      dslipmax = 1,
                      repli = 1:ncol(data$counts),
                      filter = 0,
                      param = list(m = m,
                                   sratio = sratio,
                                   contamination_rate = contamination_rate,
                                   mu = mu,
                                   nu_plus = nu_plus,
                                   nu_minus = nu_minus)){
  lambda <-  1/contamination_rate
  ll <- NULL
  for (i in repli){
    llrepli <- lltotal(data = data, genotype = genotype, repli = i, filter = filter)
    if (!is.null(llrepli)){
      llrepli$repli <- i
      ll <- rbind(ll, llrepli)
    }
  }
  ll <- rbind(data.frame(component = "Total_sample", ll = sum(ll[ll$component == "total",]$ll), repli = NA),ll)
  ll
}