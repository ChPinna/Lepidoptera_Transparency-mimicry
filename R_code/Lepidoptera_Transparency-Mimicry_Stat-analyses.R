# Mimicry drives convergence in structural and light transmission features of transparent wings in Lepidoptera - Code for analyses #

#######################
#####   RESULTS   #####
#######################

#_____________________________________________________________________________________________________
#### Convergence among co-mimics in visual appearance of transparent patches as seen by predators ####
#_____________________________________________________________________________________________________

# -------------------------------------------------
# Data, libraries and home-made function needed 
# -------------------------------------------------
library(pavo)
# Function to get colour distances calculated by the coldist function from pavo package
getColDist<-function(datatable, vision.mdl=c("VS", "UVS", "HUM"), illum=c("ws", "fs", "el", "lg"), conformation=c("transmission", "reflexion"), qcatch=c("Qi", "fi"), relative=c(FALSE, TRUE), noise = c("neural", "quantum"), subset = c("bgid", "fond_feuille"), achro){
  vs = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/sensitivity-curves/VS_shearwater.txt', header = T)
  colnames(vs) = c('wl','u','s','m','l','q')
  uvs = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/sensitivity-curves/UVS_bluetit.txt', header = T)
  colnames(uvs) = c('wl','u','s','m','l','q')
  illfs = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/light-environments/LuxFS.txt', header = T)
  ill_fs = illfs[,2]
  illel = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/light-environments/LuxEL.txt', header = T)
  ill_el = illel[,2]
  illws = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/light-environments/LuxWS.txt', header = T)
  ill_ws = illws[,2]
  illlg = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/light-environments/LuxLG.txt', header = T)
  ill_lg = illlg[,2]
  fondfeuil = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/backgrounds/FondCGfeuilles.txt', header = T)
  fond_feuille = fondfeuil[,2]
  
  if (vision.mdl=="UVS") {
    vis<-uvs 
    nvis<-c(1, 1.9, 2.7, 2.7)
    webvis<-0.1
    vision<-"uvs" #pour LWS
    visdf = vis[,c(1:5)]
    achdf = vis[,6]
    webachr=0.2
  }
  if (vision.mdl=="VS") {
    vis<-vs
    nvis<-c(1, 0.7, 1, 1.4)
    vision<-"vs"
    webvis<-0.1 #pour LWS
    visdf = vis[,c(1:5)]
    achdf = vis[,6]
    webachr=0.2
  }
  if (vision.mdl=="HUM")  {
    vis<-hum
    nvis<-c(1,16,32)
    vision<-"hum"
    webvis<-0.018 #pour LWS
    visdf = vis[,c(1:4)]
    achdf = vis[,5]
    webachr=0.11 
  }
  
  if (conformation=="transmission") {
    data<-datatable
  }
  if (conformation=="reflexion") {
    data<-as.data.frame(matrix(NA, nrow=dim(datatable)[1], ncol=dim(datatable)[2]))
    names(data)<-names(datatable)
    for (j in 2:dim(datatable)[2]){
      for (i in 1:dim(datatable)[1]){
        data[i,j]<-datatable[i,j]*datatable[i,j]*fond_feuille[i]/10000
      }
    }
    data[,1]<-datatable$wl
  }
  
  if (illum=="ws"){ 
    illdf<-ill_ws
  }
  if (illum=="fs"){ 
    illdf<-ill_fs
  }
  if (illum=="el"){ 
    illdf<-ill_el
  }
  if (illum=="lg"){ 
    illdf<-ill_lg
  }
  if (is.null(subset) == FALSE){
    if (subset == "bgid") {
      bgid<-rep(100,401)
      data <- cbind(data, as.numeric(bgid))
    }
    if (subset == "fond_feuille") {
      data <- cbind(data, fond_feuille)
    }
  }
  
  data.vismodel<-vismodel(data, visual = visdf, achromatic = achdf, illum = illdf, trans=1, qcatch=qcatch, relative=relative)
  data.coldist<-coldist(modeldata = data.vismodel, noise = noise, subset = subset, achromatic = achro, qcatch = qcatch, n = nvis, weber = webvis, weber.achro = webachr)
  
  return(data.coldist)  
}
# Function to test for convergence of traits
getConvColdist <- function(dataframe, df.class = c("coldist", "df"), data.infos, infos.tipLabel = c("TipLabel", "Tip.Label"), tree, stat = c("mean", "median", "sum"), var, family.var = c("Gaussian", "binomial"), var2fac, factor, grp, nperm, subgroup = NULL){
  # We  calculate phylogenetic distance
  if (is.null(tree) == FALSE) {
    acov=cophenetic(tree)
  }
  var2fac.nb <- which(names(data.infos) == var2fac)
  var.nb <- which(names(dataframe) == var)
  tiplab.nb <- which(names(data.infos) == infos.tipLabel)
  subgrp.nb <- which(names(data.infos) == subgroup)
  # Create dataframe with new infos
  data.temp <- dataframe
  if (df.class == "coldist"){
    ncol <- ncol(data.temp)
    for (j in 1:nrow(data.temp)){
      data.temp$Tip.Label1[j] <- data.infos[which(data.infos$SpectrumID == dataframe$patch1[j]), tiplab.nb]
      data.temp$Tip.Label2[j] <- data.infos[which(data.infos$SpectrumID == dataframe$patch2[j]), tiplab.nb]
      data.temp$Mimic1[j] <- as.character(data.infos[which(data.infos$SpectrumID == dataframe$patch1[j]), var2fac.nb])
      data.temp$Mimic2[j] <- as.character(data.infos[which(data.infos$SpectrumID == dataframe$patch2[j]), var2fac.nb])
      if (is.null(tree) == FALSE) data.temp$dist.phylo[j] <- acov[which(row.names(acov) == data.temp$Tip.Label1[j]), which(colnames(acov) == data.temp$Tip.Label2[j])]
      if (data.temp$Mimic1[j] == data.temp$Mimic2[j]){
        if (is.null(tree) == FALSE) data.temp[j,(ncol+6)] <- grp else data.temp[j,(ncol+5)] <- grp
      }
      if (data.temp$Mimic1[j] != data.temp$Mimic2[j]){
        if (is.null(tree) == FALSE) data.temp[j,(ncol+6)] <- "N" else data.temp[j,(ncol+5)] <- "N"
      }
    }
    if (is.null(tree) == FALSE) names(data.temp)[ncol+6] <- factor else names(data.temp)[ncol+5] <- factor
    if (is.null(tree) == FALSE) data.temp[,(ncol+6)] <- as.factor(as.character(data.temp[,(ncol+6)])) else data.temp[,(ncol+5)] <- as.factor(as.character(data.temp[,(ncol+5)]))
    
    if (is.null(subgroup) == FALSE){
      for (j in 1:nrow(data.temp)){
        if (is.null(tree) == FALSE) data.temp[j,(ncol+7)] <- data.infos[which(data.infos$SpectrumID == data.temp$patch1[j]), subgrp.nb]
        else data.temp[j,(ncol+6)] <- data.infos[which(data.infos$SpectrumID == data.temp$patch1[j]), subgrp.nb]
      }
      if (is.null(tree) == FALSE) names(data.temp)[ncol+7] <- subgroup else names(data.temp)[ncol+6] <- subgroup
      if (is.null(tree) == FALSE) data.temp[,(ncol+7)] <- as.factor(as.character(data.temp[,(ncol+7)])) else data.temp[,(ncol+6)] <- as.factor(as.character(data.temp[,(ncol+6)]))
    } 
  }
  
  # For a dataframe
  if (df.class == "df"){
    df.nrow <- length(unique(dataframe[,which(names(dataframe) == infos.tipLabel)]))
    if (is.null(subgroup) == TRUE) {
      data.temp <- as.data.frame(matrix(NA, nrow = df.nrow*(df.nrow-1)/2, ncol = 10)) 
      names(data.temp) <- c("Tip.Label1", "Tip.Label2", "Mimic1", "Mimic2", paste(var,".sp1", sep = ""), paste(var, ".sp2", sep = ""), "dist.phylo", "Var.Sharing", "Bino", factor)
    } 
    if(is.null(subgroup) == FALSE) {
      data.temp <- as.data.frame(matrix(NA, nrow = df.nrow*(df.nrow-1)/2, ncol = 11)) 
      names(data.temp) <- c("Tip.Label1", "Tip.Label2", "Mimic1", "Mimic2", paste(var,".sp1", sep = ""), paste(var, ".sp2", sep = ""), "dist.phylo", "Var.Sharing", "Bino", factor, subgroup)
    }
    count = 1
    for (i in 1:(df.nrow-1)){
      for (j in ((i+1):df.nrow)){
        tip1 <- data.infos[i, tiplab.nb]
        tip2 <- data.infos[j, tiplab.nb]
        mimic1 <- dataframe[which(dataframe[,tiplab.nb] == tip1), var2fac.nb]
        mimic2 <- dataframe[which(dataframe[,tiplab.nb] == tip2), var2fac.nb]
        Var.sp1 <- as.character(dataframe[which(dataframe[,tiplab.nb] == tip1), var.nb])
        Var.sp2 <- as.character(dataframe[which(dataframe[,tiplab.nb] == tip2), var.nb])
        if (is.null(tree) == FALSE) dist.phylo <- acov[which(row.names(acov) == tip1), which(colnames(acov) == tip2)]
        data.temp$Tip.Label1[count] <- tip1
        data.temp$Tip.Label2[count] <- tip2
        data.temp$Mimic1[count] <- as.character(mimic1)
        data.temp$Mimic2[count] <- as.character(mimic2)
        data.temp[count,5] <- Var.sp1
        data.temp[count,6] <- Var.sp2
        data.temp$dist.phylo[count] <- dist.phylo
        if (mimic1 == mimic2) data.temp[count, 10] <- grp else data.temp[count, 10] <- "N"
        if (Var.sp1 == Var.sp2) data.temp$`Var.Sharing`[count] <- "Y" else data.temp$`Var.Sharing`[count] <- "N"
        if (Var.sp1 == Var.sp2) data.temp$Bino[count] <- 1 else data.temp$Bino[count] <- 0
        if(is.null(subgroup) == FALSE) data.temp[count,11] <- as.character(dataframe[which(dataframe[,tiplab.nb] == tip1), subgrp.nb])
        count = count + 1
      }
    }
    var.nb <- which(names(data.temp) == "Bino")
  }
  
  factor.colnb <- which(names(data.temp) == factor)
  
  ## Before running linear model for phylogenetic correction, we calculate stat on var on the data directly and not the residuals
  # Calculation of the observed stat
  func <- as.function(get(list(stat)[[1]]))
  if (family.var == "binomial") func <- sum 
  stat.obs <- func(x = data.temp[which(data.temp[,factor.colnb] == grp), var.nb])
  # Calculation of the null distribution without phylogenetic correction
  stat.null <- NULL
  for (k in 1:nperm){
    var.rdm <- sample(data.temp[,var.nb])
    stat.rdm <- func(x = var.rdm[which(data.temp[,factor.colnb] == grp)])
    stat.null <- c(stat.null, stat.rdm)
  }
  if (family.var == "binomial"){
    p.val <- length(which(stat.null >= stat.obs))/nperm
  }
  if (family.var == "Gaussian"){
    p.val <- length(which(stat.null <= stat.obs))/nperm
  }
  ## Calculation with phylogenetic correction
  if (is.null(tree) == FALSE){
    # We run linear model
    if (family.var == "Gaussian"){
      mdl <- lm(formula = data.temp[,var.nb] ~ dist.phylo, data = data.temp)
    }
    if (family.var == "binomial"){
      mdl <- glm(formula = Bino ~ dist.phylo, family = "binomial", data = data.temp)
    }
    # We calculate the stat for observed data
    func <- as.function(get(list(stat)[[1]]))
    stat.phy.obs <- func(x = mdl$residuals[which(data.temp[,factor.colnb] == grp)])
    # We get the null distribution of the stat
    stat.phy.null <- NULL
    for (k in 1:nperm){
      res.rdm <- sample(mdl$residuals)
      stat.phy.rdm <- func(x = res.rdm[which(data.temp[,factor.colnb] == grp)])
      stat.phy.null <- c(stat.phy.null, stat.phy.rdm)
    }
    if (family.var == "binomial"){
      p.val.phy <- length(which(stat.phy.null >= stat.phy.obs))/nperm
    }
    if (family.var == "Gaussian"){
      p.val.phy <- length(which(stat.phy.null <= stat.phy.obs))/nperm
    }
  }
  
  ## Calculation for group by group
  if (is.null(subgroup) == FALSE){
    subgrp.nb <- which(names(data.temp) == subgroup)
    N.subgrp <- length(unique(data.temp[,subgrp.nb]))
    subgrp.list <- vector("list", N.subgrp)
    for (i in 1:N.subgrp){
      subgrp <- as.character(unique(data.temp[,subgrp.nb])[i])
      names(subgrp.list)[[i]] <- subgrp
      # Value of observed stat and mean residuals
      if (family.var == "binomial") func <- sum
      stat.obs.sg <- func(x = data.temp[which(data.temp[,subgrp.nb] == subgrp & data.temp[,factor.colnb] == grp), var.nb])
      func <- as.function(get(list(stat)[[1]]))
      stat.phy.obs.sg <- func(x = mdl$residuals[which(data.temp[,subgrp.nb] == subgrp & data.temp[,factor.colnb] == grp)])
      # Creation of the null distribution by randomization
      residnull.sg = vector()
      stat.null.sg = vector()
      if (family.var == "binomial") func <- sum
      for (k in 1:nperm){
        stat.null.sg <- c(stat.null.sg, func(sample(data.temp[,var.nb])[which(data.temp[,subgrp.nb] == subgrp & data.temp[,factor.colnb] == grp)]))
      }
      func <- as.function(get(list(stat)[[1]]))
      for (k in 1:nperm){
        residnull.sg <- c(residnull.sg, func(sample(mdl$residuals)[which(data.temp[,subgrp.nb] == subgrp & data.temp[,factor.colnb] == grp)]))
      }
      if (family.var == "binomial"){
        p.val.sg <- length(which(stat.null.sg >= stat.obs.sg))/nperm
        p.val.phy.sg <- length(which(residnull.sg >= stat.phy.obs.sg))/nperm
      }
      if (family.var == "Gaussian"){
        p.val.sg <- length(which(stat.null.sg <= stat.obs.sg))/nperm
        p.val.phy.sg <- length(which(residnull.sg <= stat.phy.obs.sg))/nperm
      }
      subgrp.list[[i]] <- list("withoutphycor" = list("stat.obs" = stat.obs.sg, "stat.null" = stat.null.sg, "p-value" = p.val.sg), "withphycor" = list("stat.phy.obs" = stat.phy.obs.sg, "stat.phy.null" = residnull.sg, "p-value.phy" = p.val.phy.sg))
    }
  }
  
  if (is.null(tree) == FALSE) tot.res <- list("tot" = list("datatable" = data.temp, "withoutphycor" = list("stat.obs" = stat.obs, "stat.null" = stat.null, "p-value" = p.val), "withphycor" = list("modele" = mdl, "stat.phy.obs" = stat.phy.obs, "stat.phy.null" = stat.phy.null, "p-value.phy" = p.val.phy)))
  else tot.res <- list("tot" = list("datatable" = data.temp, "withoutphycor" = list("stat.obs" = stat.obs, "stat.null" = stat.null, "p-value" = p.val)))
  if (is.null(subgroup) == FALSE) res <- list("tot" = tot.res, "subgrp" = subgrp.list)
  else res <- list("tot" = tot.res)
  
  return(res)
}
# Function to smooth spectra
spec.lo<-function(data, colnumber, span){ ## colnumber is the column number of the column used to loess
  i.in<-seq(1, dim(data)[2], by=1)
  i.in<-i.in[-colnumber]
  data.lo<-as.data.frame(data[,colnumber])
  names(data.lo)<-c("wl")
  c=2
  progbar<-txtProgressBar(min=0, max=length(i.in), style=2)
  for (i in i.in) {
    temp.lo<-as.data.frame(loess.smooth(as.vector(data[,colnumber]), as.vector(data[,i]), span=span, family="gaussian", evaluation=length(data[,colnumber])))
    data.lo$wl<-temp.lo$x
    data.lo<-cbind(data.lo, temp.lo$y)
    colnames(data.lo)[c] <- colnames(data)[i]
    c=c+1
    setTxtProgressBar(progbar,i)
  }
  return(data.lo)
}

## Raw Data
# specall.raw = dataframe with the raw spectra obtained with the getspec function from 'pavo' package
specall.raw <- read.delim(file = "C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/Lepidoptera_Transparenct-Mimicry_SpecAll-Raw.txt", header = TRUE, sep = "\t", dec = '.', check.names = FALSE)
## We smoothed the raw spectra, with the custom-made function "spec.lo"
# We fix negative value to zero
specall.cln <- procspec(specall.raw, fixneg = "zero")
# We smooth the spectra between 500 et 700nm with a span of 0.20 with the custom-made function spec.lo
specall.lo.500.700<-spec.lo(data=specall.cln[which(specall.cln$wl==500):dim(specall.cln)[1],], colnumber = 1, span=0.20)
# We prepare the datatable to smooth the entire spectrum with a span of 0.05
specall.TOLO<-rbind(specall.cln[1:which(specall.cln$wl==499),], specall.lo.500.700)
specall.lo.surlo<-spec.lo(data=specall.TOLO, colnumber = 1, span=0.05)

## Smoothed data and info about data and specimens
# specall = dataframe with the smoothed spectra (value of transmission for each wavelenght, for each spot measured)
specall <- read.delim(file = "C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/Lepidoptera_Transparenct-Mimicry_SpecAll.txt", sep = "\t", dec = ".", check.names = FALSE)
# specinfos = dataframe with the informations about spot measured and specimen
specinfos <- read.delim(file = "C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/Lepidoptera_Transparenct-Mimicry_SpecInfos.txt", header = TRUE, sep ="\t", dec = ".")
# phylogenetic tree of the species in the study
tree <- read.nexus("C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/Documents_pour-article/arbreM2_mai2020_simple.trees")

## Transformed Data
# We create a data frame for each spot (s1 to s5)
specall.spot <- vector("list", 5)
for (i in 1:5){
  toKeep <- specinfos$SpectrumID[which(specinfos$Spot == paste("s", i, sep=""))]
  specall.spot[[i]] <- specall[,c(1, which(names(specall) %in% toKeep))]
}

# We compute colour distances for each spot, for different predator visual system (UVS vs VS) and different ambiant light (forest shade and large gap)
# For UVS vision and large gap illuminant
UVS.lg.coldist <- vector('list', 5)
# For UVS vision and forest shade illuminant
UVS.fs.coldist <- vector('list', 5)
# For VS vision and large gap illuminant
VS.lg.coldist <- vector('list', 5)
# For VS vision and forest shade illuminant
VS.fs.coldist <- vector('list', 5)
# We use the home-made function getColdDist to obtain colour distance 
for (i in 1:5){
  UVS.lg.coldist[[i]] <- getColDist(datatable = specall.spot[[i]], vision.mdl = "UVS", illum = "lg", conformation = "transmission", qcatch = "fi", relative = FALSE, noise = "neural", subset = NULL, achro = TRUE)
  names(UVS.lg.coldist)[[i]] <- paste("spot s", i, sep = "")
  UVS.fs.coldist[[i]] <- getColDist(datatable = specall.spot[[i]], vision.mdl = "UVS", illum = "fs", conformation = "transmission", qcatch = "fi", relative = FALSE, noise = "neural", subset = NULL, achro = TRUE)
  names(UVS.fs.coldist)[[i]] <- paste("spot s", i, sep = "")
  VS.lg.coldist[[i]] <- getColDist(datatable = specall.spot[[i]], vision.mdl = "VS", illum = "lg", conformation = "transmission", qcatch = "fi", relative = FALSE, noise = "neural", subset = NULL, achro = TRUE)
  names(VS.lg.coldist)[[i]] <- paste("spot s", i, sep = "")
  VS.fs.coldist[[i]] <- getColDist(datatable = specall.spot[[i]], vision.mdl = "VS", illum = "fs", conformation = "transmission", qcatch = "fi", relative = FALSE, noise = "neural", subset = NULL, achro = TRUE)
  names(VS.fs.coldist)[[i]] <- paste("spot s", i, sep = "")
}


# -------------------------------------------------
# Statistical analyses
# -------------------------------------------------

#" We use the home-made function "getConvColdist" to test for convergence of optical properties
## We store the result in a list, seperately for chromatic contrast (dS) and achromatic contrast (dL)
# For UVS vision and large gap illuminant
UVS.lg.dS.phyConv <- vector('list', 5)
UVS.lg.dL.phyConv <- vector('list', 5)
# For UVS vision and forest shade illuminant
UVS.fs.dS.phyConv <- vector('list', 5)
UVS.fs.dL.phyConv <- vector('list', 5)
# For VS vision and large gap illuminant
VS.lg.dS.phyConv <- vector('list', 5)
VS.lg.dL.phyConv <- vector('list', 5)
# For VS vision and forest shade illuminant
VS.fs.dS.phyConv <- vector('list', 5)
VS.fs.dL.phyConv <- vector('list', 5)

# We use the home-made function getColdDist to obtain colour distance 
for (i in 1:5){
  UVS.lg.dS.phyConv[[i]] <- getConvColdist(dataframe = UVS.lg.coldist[[i]], df.class = "coldist", data.infos = specinfos, infos.tipLabel = "Tip.Label", tree = tree, stat = "mean", var = "dS", family.var = "Gaussian", var2fac = "mimicry.ring", factor = "comimic", grp = "Y", nperm = 10000, subgroup = "mimicry.ring")
  names(UVS.lg.dS.phyConv)[[i]] <- paste("spot s", i, sep = "")
  UVS.lg.dL.phyConv[[i]] <- getConvColdist(dataframe = UVS.lg.coldist[[i]], df.class = "coldist", data.infos = specinfos, infos.tipLabel = "Tip.Label", tree = tree, stat = "mean", var = "dL", family.var = "Gaussian", var2fac = "mimicry.ring", factor = "comimic", grp = "Y", nperm = 10000, subgroup = "mimicry.ring")
  names(UVS.lg.dL.phyConv)[[i]] <- paste("spot s", i, sep = "")
  UVS.fs.dS.phyConv[[i]] <- getConvColdist(dataframe = UVS.fs.coldist[[i]], df.class = "coldist", data.infos = specinfos, infos.tipLabel = "Tip.Label", tree = tree, stat = "mean", var = "dS", family.var = "Gaussian", var2fac = "mimicry.ring", factor = "comimic", grp = "Y", nperm = 10000, subgroup = "mimicry.ring")
  names(UVS.fs.dS.phyConv)[[i]] <- paste("spot s", i, sep = "")
  UVS.fs.dL.phyConv[[i]] <- getConvColdist(dataframe = UVS.fs.coldist[[i]], df.class = "coldist", data.infos = specinfos, infos.tipLabel = "Tip.Label", tree = tree, stat = "mean", var = "dL", family.var = "Gaussian", var2fac = "mimicry.ring", factor = "comimic", grp = "Y", nperm = 10000, subgroup = "mimicry.ring")
  names(UVS.fs.dL.phyConv)[[i]] <- paste("spot s", i, sep = "")
  VS.lg.dS.phyConv[[i]] <- getConvColdist(dataframe = VS.lg.coldist[[i]], df.class = "coldist", data.infos = specinfos, infos.tipLabel = "Tip.Label", tree = tree, stat = "mean", var = "dS", family.var = "Gaussian", var2fac = "mimicry.ring", factor = "comimic", grp = "Y", nperm = 10000, subgroup = "mimicry.ring")
  names(VS.lg.dS.phyConv)[[i]] <- paste("spot s", i, sep = "")
  VS.lg.dL.phyConv[[i]] <- getConvColdist(dataframe = VS.lg.coldist[[i]], df.class = "coldist", data.infos = specinfos, infos.tipLabel = "Tip.Label", tree = tree, stat = "mean", var = "dL", family.var = "Gaussian", var2fac = "mimicry.ring", factor = "comimic", grp = "Y", nperm = 10000, subgroup = "mimicry.ring")
  names(VS.lg.dL.phyConv)[[i]] <- paste("spot s", i, sep = "")
  VS.fs.dS.phyConv[[i]] <- getConvColdist(dataframe = VS.fs.coldist[[i]], df.class = "coldist", data.infos = specinfos, infos.tipLabel = "Tip.Label", tree = tree, stat = "mean", var = "dS", family.var = "Gaussian", var2fac = "mimicry.ring", factor = "comimic", grp = "Y", nperm = 10000, subgroup = "mimicry.ring")
  names(VS.fs.dS.phyConv)[[i]] <- paste("spot s", i, sep = "")
  VS.fs.dL.phyConv[[i]] <- getConvColdist(dataframe = VS.fs.coldist[[i]], df.class = "coldist", data.infos = specinfos, infos.tipLabel = "Tip.Label", tree = tree, stat = "mean", var = "dL", family.var = "Gaussian", var2fac = "mimicry.ring", factor = "comimic", grp = "Y", nperm = 10000, subgroup = "mimicry.ring")
  names(VS.fs.dL.phyConv)[[i]] <- paste("spot s", i, sep = "")
}

## We create a dataframe to compile results of the statistical tests computed above
coldist.phyConv.res <- as.data.frame(matrix(NA, nrow = 0, ncol =10))
names(coldist.phyConv.res) <- c("vismdl", "illum", "spot", "contrast", "stat", "group", "statObs", "p-value", "statObs.phy", "p-value.phy")
for (data in c("UVS.lg.dL.phyConv", "UVS.lg.dS.phyConv", "UVS.fs.dL.phyConv", "UVS.fs.dS.phyConv", "VS.lg.dL.phyConv", "VS.lg.dS.phyConv", "VS.fs.dL.phyConv", "VS.fs.dS.phyConv")){
  df.temp <- get(data)
  us <- gregexpr("\\.", text = data)
  vismdl <- substr(data, 1, us[[1]][1]-1)
  illum <- substr(data, us[[1]][1]+1, us[[1]][2]-1)
  contrast <- substr(data, us[[1]][2]+1, us[[1]][3]-1)
  for (i in 1:length(df.temp)){
    spot <- names(df.temp)[[i]]
    for (j in 1:(length(df.temp[[i]]$tot)+length(df.temp[[i]]$subgrp))){
      row.temp <- as.data.frame(matrix(ncol = 10))
      names(row.temp) <- names(coldist.phyConv.res)
      row.temp$vismdl <- vismdl
      row.temp$illum <- illum
      row.temp$spot <- spot
      row.temp$contrast <- contrast
      row.temp$stat <-"mean"
      if (j == 1) {
        gp.name <- names(df.temp[[i]])[[j]]
        data.temp <- df.temp[[i]][[j]][[1]]
      }
      if (j > 1) {
        gp.name <- names(df.temp[[i]][[2]])[[j-1]]
        data.temp <- df.temp[[i]][[2]][[j-1]]
      }
      row.temp$group <- gp.name
      row.temp$statObs <- data.temp$withoutphycor$stat.obs
      row.temp$`p-value` <- data.temp$withoutphycor$`p-value`
      row.temp$statObs.phy <- data.temp$withphycor$stat.phy.obs
      row.temp$`p-value.phy` <- data.temp$withphycor$`p-value.phy`
      coldist.phyConv.res <- rbind(coldist.phyConv.res, row.temp)
    }
  }
}


#_____________________________________________________________________________________________________
#### Diversity and convergence among co-mimics of structures involved in transparency ####
#_____________________________________________________________________________________________________

# _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
### Phylogenetic signal ###
# _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

# -------------------------------------------------
# Data, libraries and home-made function needed 
# -------------------------------------------------
library(phytools)
library(caper)
library(devtools)
source_url("https://raw.githubusercontent.com/mrborges23/delta_statistic/master/code.R")
# Datatable with all the structural features and mean transmittance
DataTot <- read.delim(file = "C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/Lepidoptera_Transparenct-Mimicry_Structure-Transmission-Infos-all.txt", header = TRUE, sep = "\t", dec = ".")
# We create a variable ScaleType which is the interaction between ScaleForm and ScaleInsertion
DataTot$ScaleType <- interaction(DataTot$ScaleForm, DataTot$ScaleInsertion)
DataTot$ScaleType <- droplevels(DataTot$ScaleType)
tree <- read.nexus("C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/Documents_pour-article/arbreM2_mai2020_simple.trees")

# -------------------------------------------------
# Analyses - Calculation of phylogenetic signal 
# -------------------------------------------------

## For continuous variable, Pagel's lambda and Blomberg's K
# Scale length
ScaleLength<-DataTot$`ScaleLength.mean`
names(ScaleLength)<-DataTot$TipLabel
phylosig(tree = tree, x = ScaleLength, method = "lambda", test = TRUE)
phylosig(tree = tree, x = ScaleLength, method = "K", test = TRUE)
# Scale width
ScaleWidth<-DataTot$`ScaleWidth.mean`
names(ScaleWidth)<-DataTot$TipLabel
phylosig(tree = tree, x = ScaleWidth, method = "lambda", test = TRUE)
phylosig(tree = tree, x = ScaleWidth, method = "K", test = TRUE)
# Sccale density
ScaleDensity<-DataTot$`ScaleDensity.mean`
names(ScaleDensity)<-DataTot$TipLabel
phylosig(tree = tree, x = ScaleDensity, method = "lambda", test = TRUE)
phylosig(tree = tree, x = ScaleDensity, method = "K", test = TRUE)
# Membrane thickness
MembraneThick <- DataTot$MembraneThick.mean
names(MembraneThick)<-DataTot$TipLabel
phylosig(tree = tree, x = MembraneThick, method = "lambda", test = TRUE)
phylosig(tree = tree, x = MembraneThick, method = "K", test = TRUE)
# Nanostructure density
NanoDensity<-DataTot$`NanoDensity.mean`
names(NanoDensity)<-DataTot$TipLabel
phylosig(tree = tree, x = NanoDensity, method = "lambda", test = TRUE)
phylosig(tree = tree, x = NanoDensity, method = "K", test = TRUE)
# Mean transmittance
B2mean<-DataTot$`B2tot.mean`
names(B2mean)<-DataTot$TipLabel
phylosig(tree = tree, x = B2mean, method = "lambda", test = TRUE)
phylosig(tree = tree, x = B2mean, method = "K", test = TRUE)

## For multicategorical variables (scale type and nanostructure type), delta statistic
# Nanostructure type
nano <- DataTot$NanoType
names(nano) <- DataTot$TipLabel
delta_nano <- delta(trait = nano, tree = tree, lambda0 = 0.1, se = 0.5, sim = 10000, thin = 10, burn = 100 )
# We randomize nanostructure type along the phylogeny and we calcule delta statistic to have a random distribution of delta
nano_delta_distri <- vector()
while (length(nano_delta_distri) < 1000){
  tryCatch({
    nanobis <- sample(nano, length(nano))
    names(nanobis) <- DataTot$TipLabel
    deltabis <- delta(trait = nanobis, tree = tree, lambda0 = 0.1, se = 0.5, sim = 10000, thin = 10, burn = 100)
    nano_delta_distri <- c(nano_delta_distri, deltabis)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
# Calculation of p-value as the proportion of randomizations where the value of delta is higher than the observed value of delta
nano.pval <- length(which(nano_delta_distri > delta_nano))/1000

# ScaleType
scale <- DataTot$ScaleType
names(scale) <- DataTot$TipLabel
delta_scale <- delta(trait = scale, tree = tree, lambda0 = 0.1, se = 0.5, sim = 10000, thin = 10, burn = 100)
# We randomize scale type along the phylogeny and we calcule delta statistic to have a random distribution of delta
scale_delta_distri <- vector()
while (length(scale_delta_distri) < 1000){
  tryCatch({
    scalebis <- sample(scale, length(scale))
    names(scalebis) <- DataTot$TipLabel
    deltabis <- delta(trait = scalebis, tree = tree, lambda0 = 0.1, se = 0.5, sim = 10000, thin = 10, burn = 100)
    scale_delta_distri <- c(scale_delta_distri, deltabis)}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
# Calculation of p-value as the proportion of randomizations where the value of delta is higher than the observed value of delta
scale.pval <- length(which(scale_delta_distri > delta_scale))/1000


## For categorical variables with only two categories (scale color), D statisti
#We create comparative data to be able to use the phylo.d function from caper
data.comp <- comparative.data(phy = tree, data = DataTot, names.col = TipLabel, na.omit = FALSE)
ScaleColor.D <- phylo.d(data.comp, binvar = ScaleColour, permut = 1000)



# _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
### Convergence of structural features among co-mimics ###
# _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-

# -------------------------------------------------
# Data, libraries and home-made function needed 
# -------------------------------------------------
# Datatable with all the structural features and mean transmittance
DataTot <- read.delim(file = "C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/Lepidoptera_Transparenct-Mimicry_Structure-Transmission-Infos-all.txt", header = TRUE, sep = "\t", dec = ".")
# We create the scale type variable which is the interaction between scale form and scale insertion
DataTot$ScaleType <- interaction(DataTot$ScaleForm, DataTot$ScaleInsertion)
DataTot$ScaleType <- droplevels(DataTot$ScaleType)
# We create a new variable, structural syndrome, which is the association of scale type and nanostructure type
DataTot$StructSyndrome <- interaction(DataTot$ScaleType, DataTot$NanoType)
DataTot$StructSyndrome <- droplevels(DataTot$StructSyndrome)
# custom-made function to test for convergence 
getConvColdist <- function(dataframe, df.class = c("coldist", "df"), data.infos, infos.tipLabel = c("TipLabel", "Tip.Label"), tree, stat = c("mean", "median", "sum"), var, family.var = c("Gaussian", "binomial"), var2fac, factor, grp, nperm, subgroup = NULL){
  # We  calculate phylogenetic distance
  if (is.null(tree) == FALSE) {
    acov=cophenetic(tree)
  }
  var2fac.nb <- which(names(data.infos) == var2fac)
  var.nb <- which(names(dataframe) == var)
  tiplab.nb <- which(names(data.infos) == infos.tipLabel)
  subgrp.nb <- which(names(data.infos) == subgroup)
  # Create dataframe with new infos
  data.temp <- dataframe
  if (df.class == "coldist"){
    ncol <- ncol(data.temp)
    for (j in 1:nrow(data.temp)){
      data.temp$Tip.Label1[j] <- data.infos[which(data.infos$SpectrumID == dataframe$patch1[j]), tiplab.nb]
      data.temp$Tip.Label2[j] <- data.infos[which(data.infos$SpectrumID == dataframe$patch2[j]), tiplab.nb]
      data.temp$Mimic1[j] <- as.character(data.infos[which(data.infos$SpectrumID == dataframe$patch1[j]), var2fac.nb])
      data.temp$Mimic2[j] <- as.character(data.infos[which(data.infos$SpectrumID == dataframe$patch2[j]), var2fac.nb])
      if (is.null(tree) == FALSE) data.temp$dist.phylo[j] <- acov[which(row.names(acov) == data.temp$Tip.Label1[j]), which(colnames(acov) == data.temp$Tip.Label2[j])]
      if (data.temp$Mimic1[j] == data.temp$Mimic2[j]){
        if (is.null(tree) == FALSE) data.temp[j,(ncol+6)] <- grp else data.temp[j,(ncol+5)] <- grp
      }
      if (data.temp$Mimic1[j] != data.temp$Mimic2[j]){
        if (is.null(tree) == FALSE) data.temp[j,(ncol+6)] <- "N" else data.temp[j,(ncol+5)] <- "N"
      }
    }
    if (is.null(tree) == FALSE) names(data.temp)[ncol+6] <- factor else names(data.temp)[ncol+5] <- factor
    if (is.null(tree) == FALSE) data.temp[,(ncol+6)] <- as.factor(as.character(data.temp[,(ncol+6)])) else data.temp[,(ncol+5)] <- as.factor(as.character(data.temp[,(ncol+5)]))
    
    if (is.null(subgroup) == FALSE){
      for (j in 1:nrow(data.temp)){
        if (is.null(tree) == FALSE) data.temp[j,(ncol+7)] <- data.infos[which(data.infos$SpectrumID == data.temp$patch1[j]), subgrp.nb]
        else data.temp[j,(ncol+6)] <- data.infos[which(data.infos$SpectrumID == data.temp$patch1[j]), subgrp.nb]
      }
      if (is.null(tree) == FALSE) names(data.temp)[ncol+7] <- subgroup else names(data.temp)[ncol+6] <- subgroup
      if (is.null(tree) == FALSE) data.temp[,(ncol+7)] <- as.factor(as.character(data.temp[,(ncol+7)])) else data.temp[,(ncol+6)] <- as.factor(as.character(data.temp[,(ncol+6)]))
    } 
  }
  
  # For a dataframe
  if (df.class == "df"){
    df.nrow <- length(unique(dataframe[,which(names(dataframe) == infos.tipLabel)]))
    if (is.null(subgroup) == TRUE) {
      data.temp <- as.data.frame(matrix(NA, nrow = df.nrow*(df.nrow-1)/2, ncol = 10)) 
      names(data.temp) <- c("Tip.Label1", "Tip.Label2", "Mimic1", "Mimic2", paste(var,".sp1", sep = ""), paste(var, ".sp2", sep = ""), "dist.phylo", "Var.Sharing", "Bino", factor)
    } 
    if(is.null(subgroup) == FALSE) {
      data.temp <- as.data.frame(matrix(NA, nrow = df.nrow*(df.nrow-1)/2, ncol = 11)) 
      names(data.temp) <- c("Tip.Label1", "Tip.Label2", "Mimic1", "Mimic2", paste(var,".sp1", sep = ""), paste(var, ".sp2", sep = ""), "dist.phylo", "Var.Sharing", "Bino", factor, subgroup)
    }
    count = 1
    for (i in 1:(df.nrow-1)){
      for (j in ((i+1):df.nrow)){
        tip1 <- data.infos[i, tiplab.nb]
        tip2 <- data.infos[j, tiplab.nb]
        mimic1 <- dataframe[which(dataframe[,tiplab.nb] == tip1), var2fac.nb]
        mimic2 <- dataframe[which(dataframe[,tiplab.nb] == tip2), var2fac.nb]
        Var.sp1 <- as.character(dataframe[which(dataframe[,tiplab.nb] == tip1), var.nb])
        Var.sp2 <- as.character(dataframe[which(dataframe[,tiplab.nb] == tip2), var.nb])
        if (is.null(tree) == FALSE) dist.phylo <- acov[which(row.names(acov) == tip1), which(colnames(acov) == tip2)]
        data.temp$Tip.Label1[count] <- tip1
        data.temp$Tip.Label2[count] <- tip2
        data.temp$Mimic1[count] <- as.character(mimic1)
        data.temp$Mimic2[count] <- as.character(mimic2)
        data.temp[count,5] <- Var.sp1
        data.temp[count,6] <- Var.sp2
        data.temp$dist.phylo[count] <- dist.phylo
        if (mimic1 == mimic2) data.temp[count, 10] <- grp else data.temp[count, 10] <- "N"
        if (Var.sp1 == Var.sp2) data.temp$`Var.Sharing`[count] <- "Y" else data.temp$`Var.Sharing`[count] <- "N"
        if (Var.sp1 == Var.sp2) data.temp$Bino[count] <- 1 else data.temp$Bino[count] <- 0
        if(is.null(subgroup) == FALSE) data.temp[count,11] <- as.character(dataframe[which(dataframe[,tiplab.nb] == tip1), subgrp.nb])
        count = count + 1
      }
    }
    var.nb <- which(names(data.temp) == "Bino")
  }
  
  factor.colnb <- which(names(data.temp) == factor)
  
  ## Before running linear model for phylogenetic correction, we calculate stat on var on the data directly and not the residuals
  # Calculation of the observed stat
  func <- as.function(get(list(stat)[[1]]))
  if (family.var == "binomial") func <- sum 
  stat.obs <- func(x = data.temp[which(data.temp[,factor.colnb] == grp), var.nb])
  # Calculation of the null distribution without phylogenetic correction
  stat.null <- NULL
  for (k in 1:nperm){
    var.rdm <- sample(data.temp[,var.nb])
    stat.rdm <- func(x = var.rdm[which(data.temp[,factor.colnb] == grp)])
    stat.null <- c(stat.null, stat.rdm)
  }
  if (family.var == "binomial"){
    p.val <- length(which(stat.null >= stat.obs))/nperm
  }
  if (family.var == "Gaussian"){
    p.val <- length(which(stat.null <= stat.obs))/nperm
  }
  ## Calculation with phylogenetic correction
  if (is.null(tree) == FALSE){
    # We run linear model
    if (family.var == "Gaussian"){
      mdl <- lm(formula = data.temp[,var.nb] ~ dist.phylo, data = data.temp)
    }
    if (family.var == "binomial"){
      mdl <- glm(formula = Bino ~ dist.phylo, family = "binomial", data = data.temp)
    }
    # We calculate the stat for observed data
    func <- as.function(get(list(stat)[[1]]))
    stat.phy.obs <- func(x = mdl$residuals[which(data.temp[,factor.colnb] == grp)])
    # We get the null distribution of the stat
    stat.phy.null <- NULL
    for (k in 1:nperm){
      res.rdm <- sample(mdl$residuals)
      stat.phy.rdm <- func(x = res.rdm[which(data.temp[,factor.colnb] == grp)])
      stat.phy.null <- c(stat.phy.null, stat.phy.rdm)
    }
    if (family.var == "binomial"){
      p.val.phy <- length(which(stat.phy.null >= stat.phy.obs))/nperm
    }
    if (family.var == "Gaussian"){
      p.val.phy <- length(which(stat.phy.null <= stat.phy.obs))/nperm
    }
  }
  
  ## Calculation for group by group
  if (is.null(subgroup) == FALSE){
    subgrp.nb <- which(names(data.temp) == subgroup)
    N.subgrp <- length(unique(data.temp[,subgrp.nb]))
    subgrp.list <- vector("list", N.subgrp)
    for (i in 1:N.subgrp){
      subgrp <- as.character(unique(data.temp[,subgrp.nb])[i])
      names(subgrp.list)[[i]] <- subgrp
      # Value of observed stat and mean residuals
      if (family.var == "binomial") func <- sum
      stat.obs.sg <- func(x = data.temp[which(data.temp[,subgrp.nb] == subgrp & data.temp[,factor.colnb] == grp), var.nb])
      func <- as.function(get(list(stat)[[1]]))
      stat.phy.obs.sg <- func(x = mdl$residuals[which(data.temp[,subgrp.nb] == subgrp & data.temp[,factor.colnb] == grp)])
      # Creation of the null distribution by randomization
      residnull.sg = vector()
      stat.null.sg = vector()
      if (family.var == "binomial") func <- sum
      for (k in 1:nperm){
        stat.null.sg <- c(stat.null.sg, func(sample(data.temp[,var.nb])[which(data.temp[,subgrp.nb] == subgrp & data.temp[,factor.colnb] == grp)]))
      }
      func <- as.function(get(list(stat)[[1]]))
      for (k in 1:nperm){
        residnull.sg <- c(residnull.sg, func(sample(mdl$residuals)[which(data.temp[,subgrp.nb] == subgrp & data.temp[,factor.colnb] == grp)]))
      }
      if (family.var == "binomial"){
        p.val.sg <- length(which(stat.null.sg >= stat.obs.sg))/nperm
        p.val.phy.sg <- length(which(residnull.sg >= stat.phy.obs.sg))/nperm
      }
      if (family.var == "Gaussian"){
        p.val.sg <- length(which(stat.null.sg <= stat.obs.sg))/nperm
        p.val.phy.sg <- length(which(residnull.sg <= stat.phy.obs.sg))/nperm
      }
      subgrp.list[[i]] <- list("withoutphycor" = list("stat.obs" = stat.obs.sg, "stat.null" = stat.null.sg, "p-value" = p.val.sg), "withphycor" = list("stat.phy.obs" = stat.phy.obs.sg, "stat.phy.null" = residnull.sg, "p-value.phy" = p.val.phy.sg))
    }
  }
  
  if (is.null(tree) == FALSE) tot.res <- list("tot" = list("datatable" = data.temp, "withoutphycor" = list("stat.obs" = stat.obs, "stat.null" = stat.null, "p-value" = p.val), "withphycor" = list("modele" = mdl, "stat.phy.obs" = stat.phy.obs, "stat.phy.null" = stat.phy.null, "p-value.phy" = p.val.phy)))
  else tot.res <- list("tot" = list("datatable" = data.temp, "withoutphycor" = list("stat.obs" = stat.obs, "stat.null" = stat.null, "p-value" = p.val)))
  if (is.null(subgroup) == FALSE) res <- list("tot" = tot.res, "subgrp" = subgrp.list)
  else res <- list("tot" = tot.res)
  
  return(res)
}
# phylogenetic tree with the species from the study
tree <- read.nexus("C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/Documents_pour-article/arbreM2_mai2020_simple.trees")

# -------------------------------------------------
# Analyses - Test for convergence in structures
# -------------------------------------------------
# We use the custom-made function "getConvColdist" to test for convergence of optical properties
# For scale type
ScaleType.phyConv <- getConvColdist(dataframe = DataTot, df.class = "df", data.infos = DataTot, tree = tree, stat = "mean", var = "ScaleType", family.var = "binomial", var2fac = "MimicryRing", factor = "comimic", grp = "Y", nperm = 10000, infos.tipLabel = "TipLabel", subgroup = "MimicryRing")
# For nanostructure type
NanoType.phyConv <- getConvColdist(dataframe = DataTot, df.class = "df", data.infos = DataTot, tree = tree, stat = "mean", var = "NanoType", family.var = "binomial", var2fac = "MimicryRing", factor = "comimic", grp = "Y", nperm = 10000, infos.tipLabel = "TipLabel", subgroup = "MimicryRing")
# For structural syndrome
StructTot.phyConv <- getConvColdist(dataframe = DataTot, df.class = "df", data.infos = DataTot, tree = tree, stat = "mean", var = "StructSyndrome", family.var = "binomial", var2fac = "MimicryRing", factor = "comimic", grp = "Y", nperm = 10000, infos.tipLabel = "TipLabel", subgroup = "MimicryRing")

# We create a dataframe which comprises all the results of the above convergence tests
struct.phyConv.res <- as.data.frame(matrix(NA, nrow = 0, ncol = 7))
names(struct.phyConv.res) <- c("StructuralFeature", "stat", "group", "statObs", "p-value", "statObs.phy", "p-value.phy")
for (data in c("ScaleType.phyConv", "NanoType.phyConv", "StructTot.phyConv")){
  df.temp <- get(data)
  us <- gregexpr("\\.", text = data)
  structFeat <- substr(data, 1, us[[1]][1]-1)
    for (j in 1:(length(df.temp$tot)+length(df.temp$subgrp))){
      row.temp <- as.data.frame(matrix(ncol = 7))
      names(row.temp) <- names(struct.phyConv.res)
      row.temp$StructuralFeature <- structFeat
      row.temp$stat <-"mean"
      if (j == 1) {
        gp.name <- names(df.temp)[[j]]
        data.temp <- df.temp[[j]][[1]]
      }
      if (j > 1) {
        gp.name <- names(df.temp[[2]])[[j-1]]
        data.temp <- df.temp[[2]][[j-1]]
      }
      row.temp$group <- gp.name
      row.temp$statObs <- data.temp$withoutphycor$stat.obs
      row.temp$`p-value` <- data.temp$withoutphycor$`p-value`
      row.temp$statObs.phy <- data.temp$withphycor$stat.phy.obs
      row.temp$`p-value.phy` <- data.temp$withphycor$`p-value.phy`
      struct.phyConv.res <- rbind(struct.phyConv.res, row.temp)
    }
}



#_____________________________________________________________________________________________________
#### Link between structural features and transmission properties ####
#_____________________________________________________________________________________________________

# ----------------------------
# Data and libraries  needed 
# ----------------------------
library(caper)
library(MuMIn)
DataTot <- read.delim(file = "C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/Lepidoptera_Transparenct-Mimicry_Structure-Transmission-Infos-all.txt", header = TRUE, sep = "\t", dec = ".")
tree <- read.nexus("C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/Documents_pour-article/arbreM2_mai2020_simple.trees")
# We create the ScaleType variable, which is the interaction between scale form and scale insertion
DataTot$ScaleType <- interaction(DataTot$ScaleForm, DataTot$ScaleInsertion)
DataTot$ScaleType <- droplevels(DataTot$ScaleType)
# We create the ScaleFormColour variable, which is the interaction between scale form and scale colour
DataTot$ScaleFormColour <- interaction(DataTot$ScaleForm, DataTot$ScaleColour)
DataTot$ScaleFormColour <- droplevels(DataTot$ScaleFormColour)
# We create the ScaleTypeColour variable, which is the interaction between scale form, scale insertion and scale colour
DataTot$ScaleTypeColour <- interaction(DataTot$ScaleType, DataTot$ScaleColour)
DataTot$ScaleTypeColour <- droplevels(DataTot$ScaleTypeColour)

## We create contrasts matrices for categorical variables according to what might be expected in term of transparency 
# For scale form (2 levels : piliform, lamellar)
DataTot$ScaleForm <- as.factor(as.character(DataTot$ScaleForm))
levels(DataTot$ScaleForm)
mattype<-matrix(c(-1,1),2,1) 
colnames(mattype) <- c(".P>L") # P = piliform , L = lamellar
contrasts(DataTot$ScaleForm)<- mattype
# For scale colour (2 levels : transparent or coloured)
DataTot$ScaleColour <- as.factor(as.character(DataTot$ScaleColour))
levels(DataTot$ScaleColour)
matcol<-matrix(c(-1,1),2,1) 
colnames(matcol) <- c(".T>C") # T = transparent, C = coloured
contrasts(DataTot$ScaleColour)<- matcol
# For scale insertion (2 levels : erected or flat)
DataTot$ScaleInsertion <- as.factor(as.character(DataTot$ScaleInsertion))
levels(DataTot$ScaleInsertion)
matinser<-matrix(c(-1,1),2,1) 
colnames(matinser) <- c(".E>F") # E = erected, F = flat
contrasts(DataTot$ScaleInsertion)<- matinser
# For scale colour and form (3 levels : piliform coloured, lamellar transparent and lamellar coloured)
levels(DataTot$ScaleFormColour)
mattypecol<-matrix(c(-1,2,-1,-1,0,1),3,2) 
colnames(mattypecol) <- c(".L<P", ".LC<LT") # L : lamellar , P : piliform, CL : coloured lamellar, TL : transparent lamellar
contrasts(DataTot$ScaleFormColour)<- mattypecol
# For scale type, which is the interaction between scale form and scale insertion (3 levels : piliform, lamellar erected and lamellar flat)
levels(DataTot$ScaleType)
mattypeinser<-matrix(c(-1,2,-1,1,0,-1),3,2) # 2 correspond au nombre de niveau de ma variable; 1 correspond au nombre de contraste qu'on peut faire, soit n-1=1 ici
colnames(mattypeinser) <- c(".L<P", ".LF<LE") # L : lamellar, P : piliform; LF : lamellar flat; LE : lamellar erected
contrasts(DataTot$ScaleType)<- mattypeinser
# For scale type and color, which is the interaction between scale type and scale colour (5 levels : piliform, lamellar flat coloured, lamellar flat transparent, lamellar erected coloured, lamellar erected transparent)
levels(DataTot$ScaleTypeColour)
mattypecolinser<-matrix(c(-1,4,-1,-1,-1,1,0,-1,1,-1,0,0,-1,0,1,-1,0,0,1,0),5,4)
colnames(mattypecolinser) <- c(".L<P", ".LF<LE", ".LFC<LFT", ".LEC<LET") # L : lamellar, P : piliform, LF : lamellar flat; LE : lamellar erected;  LFC : lamellar flat coloured; LFT : lamellar flat transparent; LEC : lamellar erected coloured; LET : lamellar erected transparent
contrasts(DataTot$ScaleTypeColour)<- mattypecolinser

# We create a comparative.data object to apply caper::pgls function afterwards
data.comp <- comparative.data(phy = tree, data = DataTot, vcv = TRUE, names.col = TipLabel, na.omit = FALSE)


# -------------------------------------------------
# Analyses - pgls models
# -------------------------------------------------

## We use the dredge function to test all models possible from a complete models
# We define the complete model as a pgls model
mdl.tot.pgls<-pgls(formula = B2.mean ~ 1 + ScaleType + ScaleColour + ScaleDensity.mean + ScaleLength.mean + ScaleWidth.mean + MembraneThick.mean + NanoType + NanoDensity.mean + ScaleType:ScaleDensity.mean + ScaleDensity.mean:NanoDensity.mean + ScaleDensity.mean:ScaleLength.mean:ScaleWidth.mean, data=data.comp, lambda = "ML")
# On lance la fonction MuMIn::dredge qui permet de tester de manière exhaustive tous les modèles en empêchant certaines coexistence de variables
mdl.tot.pgls.AICc.subset<-dredge(mdl.tot.pgls, beta="none" , evaluate = TRUE, rank="AICc", subset = ! c(ScaleDensity.mean && NanoType, ScaleLength.mean && NanoType, ScaleLength.mean && ScaleType, ScaleLength.mean && ScaleColour, ScaleWidth.mean && ScaleType, ScaleWidth.mean && ScaleColour, ScaleWidth.mean && NanoType, MembraneThick.mean && NanoType, NanoDensity.mean && ScaleType, NanoDensity.mean && NanoType))
# On ne considère que les meilleurs modèles (d'après leur AICc) qui ont un delta(AICc)<2 avec la fonction MuMIn::get.models
mdl.selec.pgls.AICc.subset<-get.models(mdl.tot.pgls.AICc.subset, subset = delta < 2)



######################################
#####   SUPPLEMENTARY MATERIAL   #####
######################################

#_____________________________________________________________________________________________________
#### Spectrophotometry - Repeatability ####
#_____________________________________________________________________________________________________

# -------------------------------------------------
# Data, libraries and home-made function needed 
# -------------------------------------------------
library(pavo)
library(rptR)
# Function to get colour distances calculated by the coldist function from pavo package
getColDist<-function(datatable, vision.mdl=c("VS", "UVS", "HUM"), illum=c("ws", "fs", "el", "lg"), conformation=c("transmission", "reflexion"), qcatch=c("Qi", "fi"), relative=c(FALSE, TRUE), noise = c("neural", "quantum"), subset = c("bgid", "fond_feuille"), achro, lim=c(300, 700)){
  vs = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/sensitivity-curves/VS_shearwater.txt', header = T)
  colnames(vs) = c('wl','u','s','m','l','q')
  uvs = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/sensitivity-curves/UVS_bluetit.txt', header = T)
  colnames(uvs) = c('wl','u','s','m','l','q')
  illfs = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/light-environments/LuxFS.txt', header = T)
  ill_fs = illfs[which(illfs[,1] >= lim[1] & illfs[,1] <= lim[2]),2]
  illel = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/light-environments/LuxEL.txt', header = T)
  ill_el = illel[which(illel[,1] >= lim[1] & illel[,1] <= lim[2]),2]
  illws = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/light-environments/LuxWS.txt', header = T)
  ill_ws = illws[which(illws[,1] >= lim[1] & illws[,1] <= lim[2]),2]
  illlg = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/light-environments/LuxLG.txt', header = T)
  ill_lg = illlg[which(illlg[,1] >= lim[1] & illlg[,1] <= lim[2]),2]
  fondfeuil = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/backgrounds/FondCGfeuilles.txt', header = T)
  fond_feuille = fondfeuil[which(fondfeuil[,1] >= lim[1] & fondfeuil[,1] <= lim[2]),2]
  
  if (vision.mdl=="UVS") {
    vis<-uvs[which(uvs[,1] >= lim[1] & uvs[,1] <= lim[2]),]
    nvis<-c(1, 1.9, 2.7, 2.7)
    webvis<-0.1
    vision<-"uvs" #pour LWS
    visdf = vis[,c(1:5)]
    achdf = vis[,6]
    webachr=0.2
  }
  if (vision.mdl=="VS") {
    vis<-vs[which(vs[,1] >= lim[1] & vs[,1] <= lim[2]),]
    nvis<-c(1, 0.7, 1, 1.4)
    vision<-"vs"
    webvis<-0.1 #pour LWS
    visdf = vis[,c(1:5)]
    achdf = vis[,6]
    webachr=0.2
  }
  if (vision.mdl=="HUM")  {
    vis<-hum[which(hum[,1] >= lim[1] & hum[,1] <= lim[2]),]
    nvis<-c(1,16,32)
    vision<-"hum"
    webvis<-0.018 #pour LWS
    visdf = vis[,c(1:4)]
    achdf = vis[,5]
    webachr=0.11 
  }
  
  if (conformation=="transmission") {
    data<-datatable
  }
  if (conformation=="reflexion") {
    data<-as.data.frame(matrix(NA, nrow=dim(datatable)[1], ncol=dim(datatable)[2]))
    names(data)<-names(datatable)
    for (j in 2:dim(datatable)[2]){
      for (i in 1:dim(datatable)[1]){
        data[i,j]<-datatable[i,j]*datatable[i,j]*fond_feuille[i]/10000
      }
    }
    data[,1]<-datatable$wl
  }
  
  if (illum=="ws"){ 
    illdf<-ill_ws
  }
  if (illum=="fs"){ 
    illdf<-ill_fs
  }
  if (illum=="el"){ 
    illdf<-ill_el
  }
  if (illum=="lg"){ 
    illdf<-ill_lg
  }
  if (is.null(subset) == FALSE){
    if (subset == "bgid") {
      bgid<-rep(100,(lim[2]-lim[1]+1))
      data <- cbind(data, as.numeric(bgid))
    }
    if (subset == "fond_feuille") {
      data <- cbind(data, fond_feuille)
    }
  }
  
  data.vismodel<-vismodel(data, visual = visdf, achromatic = achdf, illum = illdf, trans=1, qcatch=qcatch, relative=relative)
  data.coldist<-coldist(modeldata = data.vismodel, noise = noise, subset = subset, achromatic = achro, qcatch = qcatch, n = nvis, weber = webvis, weber.achro = webachr)
  
  return(data.coldist)  
}
# Function to test for convergence of traits
getConvColdist <- function(dataframe, df.class = c("coldist", "df"), data.infos, infos.tipLabel = c("TipLabel", "Tip.Label"), tree, stat = c("mean", "median", "sum"), var, family.var = c("Gaussian", "binomial"), var2fac, factor, grp, nperm, subgroup = NULL){
  # We  calculate phylogenetic distance
  if (is.null(tree) == FALSE) {
    acov=cophenetic(tree)
  }
  var2fac.nb <- which(names(data.infos) == var2fac)
  var.nb <- which(names(dataframe) == var)
  tiplab.nb <- which(names(data.infos) == infos.tipLabel)
  subgrp.nb <- which(names(data.infos) == subgroup)
  # Create dataframe with new infos
  data.temp <- dataframe
  if (df.class == "coldist"){
    ncol <- ncol(data.temp)
    for (j in 1:nrow(data.temp)){
      data.temp$Tip.Label1[j] <- data.infos[which(data.infos$SpectrumID == dataframe$patch1[j]), tiplab.nb]
      data.temp$Tip.Label2[j] <- data.infos[which(data.infos$SpectrumID == dataframe$patch2[j]), tiplab.nb]
      data.temp$Mimic1[j] <- as.character(data.infos[which(data.infos$SpectrumID == dataframe$patch1[j]), var2fac.nb])
      data.temp$Mimic2[j] <- as.character(data.infos[which(data.infos$SpectrumID == dataframe$patch2[j]), var2fac.nb])
      if (is.null(tree) == FALSE) data.temp$dist.phylo[j] <- acov[which(row.names(acov) == data.temp$Tip.Label1[j]), which(colnames(acov) == data.temp$Tip.Label2[j])]
      if (data.temp$Mimic1[j] == data.temp$Mimic2[j]){
        if (is.null(tree) == FALSE) data.temp[j,(ncol+6)] <- grp else data.temp[j,(ncol+5)] <- grp
      }
      if (data.temp$Mimic1[j] != data.temp$Mimic2[j]){
        if (is.null(tree) == FALSE) data.temp[j,(ncol+6)] <- "N" else data.temp[j,(ncol+5)] <- "N"
      }
    }
    if (is.null(tree) == FALSE) names(data.temp)[ncol+6] <- factor else names(data.temp)[ncol+5] <- factor
    if (is.null(tree) == FALSE) data.temp[,(ncol+6)] <- as.factor(as.character(data.temp[,(ncol+6)])) else data.temp[,(ncol+5)] <- as.factor(as.character(data.temp[,(ncol+5)]))
    
    if (is.null(subgroup) == FALSE){
      for (j in 1:nrow(data.temp)){
        if (is.null(tree) == FALSE) data.temp[j,(ncol+7)] <- data.infos[which(data.infos$SpectrumID == data.temp$patch1[j]), subgrp.nb]
        else data.temp[j,(ncol+6)] <- data.infos[which(data.infos$SpectrumID == data.temp$patch1[j]), subgrp.nb]
      }
      if (is.null(tree) == FALSE) names(data.temp)[ncol+7] <- subgroup else names(data.temp)[ncol+6] <- subgroup
      if (is.null(tree) == FALSE) data.temp[,(ncol+7)] <- as.factor(as.character(data.temp[,(ncol+7)])) else data.temp[,(ncol+6)] <- as.factor(as.character(data.temp[,(ncol+6)]))
    } 
  }
  
  # For a dataframe
  if (df.class == "df"){
    df.nrow <- length(unique(dataframe[,which(names(dataframe) == infos.tipLabel)]))
    if (is.null(subgroup) == TRUE) {
      data.temp <- as.data.frame(matrix(NA, nrow = df.nrow*(df.nrow-1)/2, ncol = 10)) 
      names(data.temp) <- c("Tip.Label1", "Tip.Label2", "Mimic1", "Mimic2", paste(var,".sp1", sep = ""), paste(var, ".sp2", sep = ""), "dist.phylo", "Var.Sharing", "Bino", factor)
    } 
    if(is.null(subgroup) == FALSE) {
      data.temp <- as.data.frame(matrix(NA, nrow = df.nrow*(df.nrow-1)/2, ncol = 11)) 
      names(data.temp) <- c("Tip.Label1", "Tip.Label2", "Mimic1", "Mimic2", paste(var,".sp1", sep = ""), paste(var, ".sp2", sep = ""), "dist.phylo", "Var.Sharing", "Bino", factor, subgroup)
    }
    count = 1
    for (i in 1:(df.nrow-1)){
      for (j in ((i+1):df.nrow)){
        tip1 <- data.infos[i, tiplab.nb]
        tip2 <- data.infos[j, tiplab.nb]
        mimic1 <- dataframe[which(dataframe[,tiplab.nb] == tip1), var2fac.nb]
        mimic2 <- dataframe[which(dataframe[,tiplab.nb] == tip2), var2fac.nb]
        Var.sp1 <- as.character(dataframe[which(dataframe[,tiplab.nb] == tip1), var.nb])
        Var.sp2 <- as.character(dataframe[which(dataframe[,tiplab.nb] == tip2), var.nb])
        if (is.null(tree) == FALSE) dist.phylo <- acov[which(row.names(acov) == tip1), which(colnames(acov) == tip2)]
        data.temp$Tip.Label1[count] <- tip1
        data.temp$Tip.Label2[count] <- tip2
        data.temp$Mimic1[count] <- as.character(mimic1)
        data.temp$Mimic2[count] <- as.character(mimic2)
        data.temp[count,5] <- Var.sp1
        data.temp[count,6] <- Var.sp2
        data.temp$dist.phylo[count] <- dist.phylo
        if (mimic1 == mimic2) data.temp[count, 10] <- grp else data.temp[count, 10] <- "N"
        if (Var.sp1 == Var.sp2) data.temp$`Var.Sharing`[count] <- "Y" else data.temp$`Var.Sharing`[count] <- "N"
        if (Var.sp1 == Var.sp2) data.temp$Bino[count] <- 1 else data.temp$Bino[count] <- 0
        if(is.null(subgroup) == FALSE) data.temp[count,11] <- as.character(dataframe[which(dataframe[,tiplab.nb] == tip1), subgrp.nb])
        count = count + 1
      }
    }
    var.nb <- which(names(data.temp) == "Bino")
  }
  
  factor.colnb <- which(names(data.temp) == factor)
  
  ## Before running linear model for phylogenetic correction, we calculate stat on var on the data directly and not the residuals
  # Calculation of the observed stat
  func <- as.function(get(list(stat)[[1]]))
  if (family.var == "binomial") func <- sum 
  stat.obs <- func(x = data.temp[which(data.temp[,factor.colnb] == grp), var.nb])
  # Calculation of the null distribution without phylogenetic correction
  stat.null <- NULL
  for (k in 1:nperm){
    var.rdm <- sample(data.temp[,var.nb])
    stat.rdm <- func(x = var.rdm[which(data.temp[,factor.colnb] == grp)])
    stat.null <- c(stat.null, stat.rdm)
  }
  if (family.var == "binomial"){
    p.val <- length(which(stat.null >= stat.obs))/nperm
  }
  if (family.var == "Gaussian"){
    p.val <- length(which(stat.null <= stat.obs))/nperm
  }
  ## Calculation with phylogenetic correction
  if (is.null(tree) == FALSE){
    # We run linear model
    if (family.var == "Gaussian"){
      mdl <- lm(formula = data.temp[,var.nb] ~ dist.phylo, data = data.temp)
    }
    if (family.var == "binomial"){
      mdl <- glm(formula = Bino ~ dist.phylo, family = "binomial", data = data.temp)
    }
    # We calculate the stat for observed data
    func <- as.function(get(list(stat)[[1]]))
    stat.phy.obs <- func(x = mdl$residuals[which(data.temp[,factor.colnb] == grp)])
    # We get the null distribution of the stat
    stat.phy.null <- NULL
    for (k in 1:nperm){
      res.rdm <- sample(mdl$residuals)
      stat.phy.rdm <- func(x = res.rdm[which(data.temp[,factor.colnb] == grp)])
      stat.phy.null <- c(stat.phy.null, stat.phy.rdm)
    }
    if (family.var == "binomial"){
      p.val.phy <- length(which(stat.phy.null >= stat.phy.obs))/nperm
    }
    if (family.var == "Gaussian"){
      p.val.phy <- length(which(stat.phy.null <= stat.phy.obs))/nperm
    }
  }
  
  ## Calculation for group by group
  if (is.null(subgroup) == FALSE){
    subgrp.nb <- which(names(data.temp) == subgroup)
    N.subgrp <- length(unique(data.temp[,subgrp.nb]))
    subgrp.list <- vector("list", N.subgrp)
    for (i in 1:N.subgrp){
      subgrp <- as.character(unique(data.temp[,subgrp.nb])[i])
      names(subgrp.list)[[i]] <- subgrp
      # Value of observed stat and mean residuals
      if (family.var == "binomial") func <- sum
      stat.obs.sg <- func(x = data.temp[which(data.temp[,subgrp.nb] == subgrp & data.temp[,factor.colnb] == grp), var.nb])
      func <- as.function(get(list(stat)[[1]]))
      stat.phy.obs.sg <- func(x = mdl$residuals[which(data.temp[,subgrp.nb] == subgrp & data.temp[,factor.colnb] == grp)])
      # Creation of the null distribution by randomization
      residnull.sg = vector()
      stat.null.sg = vector()
      if (family.var == "binomial") func <- sum
      for (k in 1:nperm){
        stat.null.sg <- c(stat.null.sg, func(sample(data.temp[,var.nb])[which(data.temp[,subgrp.nb] == subgrp & data.temp[,factor.colnb] == grp)]))
      }
      func <- as.function(get(list(stat)[[1]]))
      for (k in 1:nperm){
        residnull.sg <- c(residnull.sg, func(sample(mdl$residuals)[which(data.temp[,subgrp.nb] == subgrp & data.temp[,factor.colnb] == grp)]))
      }
      if (family.var == "binomial"){
        p.val.sg <- length(which(stat.null.sg >= stat.obs.sg))/nperm
        p.val.phy.sg <- length(which(residnull.sg >= stat.phy.obs.sg))/nperm
      }
      if (family.var == "Gaussian"){
        p.val.sg <- length(which(stat.null.sg <= stat.obs.sg))/nperm
        p.val.phy.sg <- length(which(residnull.sg <= stat.phy.obs.sg))/nperm
      }
      subgrp.list[[i]] <- list("withoutphycor" = list("stat.obs" = stat.obs.sg, "stat.null" = stat.null.sg, "p-value" = p.val.sg), "withphycor" = list("stat.phy.obs" = stat.phy.obs.sg, "stat.phy.null" = residnull.sg, "p-value.phy" = p.val.phy.sg))
    }
  }
  
  if (is.null(tree) == FALSE) tot.res <- list("tot" = list("datatable" = data.temp, "withoutphycor" = list("stat.obs" = stat.obs, "stat.null" = stat.null, "p-value" = p.val), "withphycor" = list("modele" = mdl, "stat.phy.obs" = stat.phy.obs, "stat.phy.null" = stat.phy.null, "p-value.phy" = p.val.phy)))
  else tot.res <- list("tot" = list("datatable" = data.temp, "withoutphycor" = list("stat.obs" = stat.obs, "stat.null" = stat.null, "p-value" = p.val)))
  if (is.null(subgroup) == FALSE) res <- list("tot" = tot.res, "subgrp" = subgrp.list)
  else res <- list("tot" = tot.res)
  
  return(res)
}
# Custom-made function to smooth spectra
spec.lo<-function(data, colnumber, span){ ## colnumber is the column number of the column used to loess
  i.in<-seq(1, dim(data)[2], by=1)
  i.in<-i.in[-colnumber]
  data.lo<-as.data.frame(data[,colnumber])
  names(data.lo)<-c("wl")
  c=2
  progbar<-txtProgressBar(min=0, max=length(i.in), style=2)
  for (i in i.in) {
    temp.lo<-as.data.frame(loess.smooth(as.vector(data[,colnumber]), as.vector(data[,i]), span=span, family="gaussian", evaluation=length(data[,colnumber])))
    data.lo$wl<-temp.lo$x
    data.lo<-cbind(data.lo, temp.lo$y)
    colnames(data.lo)[c] <- colnames(data)[i]
    c=c+1
    setTxtProgressBar(progbar,i)
  }
  return(data.lo)
}

## Raw data
spec4Repet.raw <- read.delim(file = "C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/Lepidoptera_Transparenct-Mimicry_Spec-4Repet-Raw.txt", header = TRUE, sep = "\t", dec = ".", check.names = FALSE)
# We smooth data, transform it in a rspect object and correct for negative value
spec4Repet <- spec.lo(data = spec4Repet.raw, colnumber = 1, span = 0.02)
spec4Repet <- as.rspec(spec4Repet)
spec4Repet <- procspec(rspecdata = spec4Repet, fixneg = "zero")
# Datatable with information about spectra and specimens
infos4Repet <- read.delim(file = "C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/Lepidoptera_Transparenct-Mimicry_SpecInfos-4Repet.txt", header = TRUE, sep = "\t", check.names = FALSE)
## We extract the interesting shape descriptors (brightness B2 and chroma S8) from smoothed spectra
shape4Repet <- summary(spec4Repet, subset = c("B2", "S8"))
shape4Repet$SpectrumID <- row.names(shape4Repet)
optics4Repet <- merge(shape4Repet, infos4Repet, by = "SpectrumID")

# --------------------------------------------------------
# Biological repeatability - intraspecific repeatability
# --------------------------------------------------------
### First we consider only the first measurement for each spot of each specimen for each species
optics4Repet.bio <- subset(optics4Repet, Repet == "R01")
# We create a variable SpeSpot which is the interaction of species and spot on the forewing
optics4Repet.bio$SpeSpot <- interaction(optics4Repet.bio$SpeShort, optics4Repet.bio$Veine)

### Calculation of the value of R with rpt function
# We look at the distribution of B2
hist(optics4Repet.bio$B2, nclass=30)
hist(optics4Repet.bio$S8, nclass=30)
# Even if the distribution is not Gaussian, we assumre that this is the closest datatype for our kind of data
# For mean transmission over 400-700 (B2)
rep.bio.B2 <- rpt(formula = B2 ~ (1|SpeShort) + (1|SpeSpot), grname = c("SpeShort", "SpeSpot"), data = optics4Repet.bio, datatype = "Gaussian", nboot = 1000, nperm = 1000)
summary(rep.bio.B2)
# For chroma (S8)
rep.bio.S8 <- rpt(formula = S8 ~ (1|SpeShort) + (1|SpeSpot), grname = c("SpeShort", "SpeSpot"), data = optics4Repet.bio, datatype = "Gaussian", nboot = 1000, nperm= 1000)
summary(rep.bio.S8)


### Calculation of median coefficient of variation for SpeSpot
# For mean transmission over 400-700nm (B2)
B2.bio.CV <- as.data.frame(matrix(NA, nrow = 19*5, ncol = 5 ))
names(B2.bio.CV) <- c("SpeSpot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(B2.bio.CV)){
  spespot <- levels(as.factor(as.character(optics4Repet.bio$SpeSpot)))[i]
  MEAN <- mean(optics4Repet.bio$B2[which(optics4Repet.bio$SpeSpot == spespot)])
  SD <- sd(optics4Repet.bio$B2[which(optics4Repet.bio$SpeSpot == spespot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  B2.bio.CV$SpeSpot[i] <- spespot
  B2.bio.CV$mean[i] <- MEAN
  B2.bio.CV$sd[i] <- SD
  B2.bio.CV$CV[i] <- CV
  B2.bio.CV$RSD[i] <- RSD*100
}
median(B2.bio.CV$RSD)

# For chroma (S8)
S8.bio.CV <- as.data.frame(matrix(NA, nrow = 19*5, ncol = 5 ))
names(S8.bio.CV) <- c("SpeSpot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(S8.bio.CV)){
  spespot <- levels(as.factor(as.character(optics4Repet.bio$SpeSpot)))[i]
  MEAN <- mean(optics4Repet.bio$S8[which(optics4Repet.bio$SpeSpot == spespot)])
  SD <- sd(optics4Repet.bio$S8[which(optics4Repet.bio$SpeSpot == spespot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  S8.bio.CV$SpeSpot[i] <- spespot
  S8.bio.CV$mean[i] <- MEAN
  S8.bio.CV$sd[i] <- SD
  S8.bio.CV$CV[i] <- CV
  S8.bio.CV$RSD[i] <- RSD*100
}
median(S8.bio.CV$RSD)


### Test of the similarity between conspecific individuals
# We create one dataframe by spot (on the FW) considered
optics4Repet.list <- vector("list", 5)
for (i in 1:5){
  toKeep <- optics4Repet$SpectrumID[which(optics4Repet$Repet == "R01" & optics4Repet$Veine == paste("v", i, sep=""))]
  optics4Repet.list[[i]] <- spec4Repet[,c(1, which(names(spec4Repet) %in% toKeep))]
}
# We calculate colour distances with the custom-made function getColdDist
optics4Repet.coldist <- vector("list", 5)
for (i in 1:5){
  optics4Repet.coldist[[i]] <- getColDist(datatable = optics4Repet.list[[i]], vision.mdl = "UVS", illum = "lg", conformation = "transmission", qcatch = "fi", relative = FALSE, noise = "neural", subset = NULL, achro = TRUE, lim = c(400,700))
}

# We now calculate with the custom-made function the similarity
optics4Repet.dS <- vector("list", 5)
optics4Repet.dL <- vector("list", 5)
for (i in 1:5){
  optics4Repet.dS[[i]] <- getConvColdist(dataframe = optics4Repet.coldist[[i]], df.class = "coldist", data.infos = optics4Repet.bio, infos.tipLabel = "SpectrumID", tree = NULL, stat = "mean", var = "dS", family.var = "Gaussian", var2fac = "SpeShort", factor = "conspecific", grp = "Y", nperm = 10000, subgroup = NULL)
  optics4Repet.dL[[i]] <- getConvColdist(dataframe = optics4Repet.coldist[[i]], df.class = "coldist", data.infos = optics4Repet.bio, infos.tipLabel = "SpectrumID", tree = NULL, stat = "mean", var = "dL", family.var = "Gaussian", var2fac = "SpeShort", factor = "conspecific", grp = "Y", nperm = 10000, subgroup = NULL)
}

# We create a dataframe with the results of the similarity test
optics4Repet.res <- as.data.frame(matrix(NA, nrow = 0, ncol =8))
names(optics4Repet.res) <- c("vismdl", "illum", "spot", "contrast", "stat", "group", "statObs", "p-value")
for (data in c("optics4Repet.dS", "optics4Repet.dL")){
  df.temp <- get(data)
  us <- gregexpr("\\.", text = data)
  vismdl <- "UVS"
  illum <- "lg"
  contrast <- substr(data, us[[1]][1]+1, nchar(data))
  for (i in 1:length(df.temp)){
    spot <- paste("s", i, sep = "")
    for (j in 1:(length(df.temp[[i]]$tot))){
      row.temp <- as.data.frame(matrix(ncol = 8))
      names(row.temp) <- names(optics4Repet.res)
      row.temp$vismdl <- vismdl
      row.temp$illum <- illum
      row.temp$spot <- spot
      row.temp$contrast <- contrast
      row.temp$stat <-"mean"
      if (j == 1) {
        gp.name <- names(df.temp[[i]])[[j]]
        data.temp <- df.temp[[i]][[j]][[1]]
      }
      if (j > 1) {
        gp.name <- names(df.temp[[i]][[2]])[[j-1]]
        data.temp <- df.temp[[i]][[2]][[j-1]]
      }
      row.temp$group <- gp.name
      row.temp$statObs <- data.temp$withoutphycor$stat.obs
      row.temp$`p-value` <- data.temp$withoutphycor$`p-value`
      optics4Repet.res <- rbind(optics4Repet.res, row.temp)
    }
  }
}


# --------------------------------------------------------
# Technical repeatability - repeatability of the measure
# --------------------------------------------------------
### First we consider only the spots which has been measure 3 times
toKeep <- unique(optics4Repet$SpeShort[which(optics4Repet$Repet == "R03")])
optics4Repet.tech <- optics4Repet[which(optics4Repet$SpeShort %in% toKeep),]
# We create a variable Spot which is the interaction of individual and spot on the forewing
optics4Repet.tech$Spot <- interaction(optics4Repet.tech$SpecimenID, optics4Repet.tech$Veine)

### Calculation of the value of R with rpt function
# We look at the distribution of B2 and S8
hist(optics4Repet.tech$B2, nclass=30)
hist(optics4Repet.tech$S8, nclass=30)
# Even if the distribution is not Gaussian, we assumre that this is the closest datatype for our kind of data
# For mean transmission over 400-700 (B2)
rep.tech.B2 <- rpt(formula = B2 ~ (1|SpeShort) + (1|SpecimenID) + (1|Spot), grname = c("SpeShort", "SpecimenID", "Spot"), data = optics4Repet.tech, datatype = "Gaussian", nboot = 1000, nperm = 1000)
summary(rep.tech.B2)
# For chroma (S8)
rep.tech.S8 <- rpt(formula = S8 ~ (1|SpeShort) + (1|SpecimenID) + (1|Spot), grname = c("SpeShort", "SpecimenID", "Spot"), data = optics4Repet.tech, datatype = "Gaussian", nboot = 1000, nperm= 1000)
summary(rep.tech.S8)


### Calculation of median coefficient of variation for SpeSpot
# For mean transmission over 400-700nm (B2)
B2.tech.CV <- as.data.frame(matrix(NA, nrow = 32*5, ncol = 5 ))
names(B2.tech.CV) <- c("Spot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(B2.tech.CV)){
  spot <- levels(as.factor(as.character(optics4Repet.tech$Spot)))[i]
  MEAN <- mean(optics4Repet.tech$B2[which(optics4Repet.tech$Spot == spot)])
  SD <- sd(optics4Repet.tech$B2[which(optics4Repet.tech$Spot == spot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  B2.tech.CV$Spot[i] <- spot
  B2.tech.CV$mean[i] <- MEAN
  B2.tech.CV$sd[i] <- SD
  B2.tech.CV$CV[i] <- CV
  B2.tech.CV$RSD[i] <- RSD*100
}
median(B2.tech.CV$RSD, na.rm = TRUE)

# For chroma (S8)
S8.tech.CV <- as.data.frame(matrix(NA, nrow = 32*5, ncol = 5 ))
names(S8.tech.CV) <- c("Spot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(S8.tech.CV)){
  spot <- levels(as.factor(as.character(optics4Repet.tech$Spot)))[i]
  MEAN <- mean(optics4Repet.tech$S8[which(optics4Repet.tech$Spot == spot)])
  SD <- sd(optics4Repet.tech$S8[which(optics4Repet.tech$Spot == spot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  S8.tech.CV$SpeSpot[i] <- spespot
  S8.tech.CV$mean[i] <- MEAN
  S8.tech.CV$sd[i] <- SD
  S8.tech.CV$CV[i] <- CV
  S8.tech.CV$RSD[i] <- RSD*100
}
median(S8.tech.CV$RSD, na.rm= TRUE)


#_____________________________________________________________________________________________________
#### Structural features - Repeatability ####
#_____________________________________________________________________________________________________

# -------------------------------------------------
# Data, libraries and home-made function needed 
# -------------------------------------------------
library(rptR)
## Data
Density<-read.delim(file = "C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/var-intra_densite_final.txt", sep="\t", header=T)
Scales<-read.delim("C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/variation-intra_longueur-largeur.txt", sep="\t", header=T)
ScalesBi<-read.delim("C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/var-intra_poils_bifides.txt", sep="\t", header=T)
ScalesMono<-read.delim("C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/var-intra_poils-monofides.txt", sep="\t", header=T)

# --------------------------------------------------------
# Biological repeatability - intraspecific repeatability
# --------------------------------------------------------
# We create a variable SpeSpot which is the interaction of species and spot on the forewing
Density$SpeSpot <- interaction(Density$Genus, Density$Zone)
ScalesBi$SpeSpot <- interaction(ScalesBi$Genus, ScalesBi$Veine)
ScalesMono$SpeSpot <- interaction(ScalesMono$Genus, ScalesMono$Veine)

### Calculation of the value of R with rpt function
# We look at the distribution of scale density, scale length and scale width
hist(Density$densite.pedicelle, nclass = 30)
hist(ScalesBi$LongueurTot1, nclass = 30)
hist(ScalesBi$Largeur.Base, nclass = 30)
hist(ScalesMono$LongueurTot1, nclass = 30)
hist(ScalesMono$Largeur.Base, nclass = 30)
# We consider that the variable are normally distributed because it is at least the closest distribution available in rptR
## For scale density
rep.bio.density <- rpt(formula = densite.pedicelle ~ (1|SpeSpot), grname = c("SpeSpot"), data = Density, datatype = "Gaussian", nboot = 1000, nperm = 1000)
summary(rep.bio.density)
## For scale length
# For bifid piliform scale
rep.bio.length.bi <- rpt(formula = LongueurTot1 ~ (1|SpeSpot), grname = c("SpeSpot"), data = ScalesBi, datatype = "Gaussian", nboot = 1000, nperm= 1000)
summary(rep.bio.length.bi)
# For monofid piliform scale
rep.bio.length.mono <- rpt(formula = LongueurTot1 ~ (1|SpeSpot), grname = c("SpeSpot"), data = ScalesMono, datatype = "Gaussian", nboot = 1000, nperm= 1000)
summary(rep.bio.length.mono)
## For scale width
# For bifid piliform scale
rep.bio.width.bi <- rpt(formula = Largeur.Base ~ (1|SpeSpot), grname = c("SpeSpot"), data = ScalesBi, datatype = "Gaussian", nboot = 1000, nperm= 1000)
summary(rep.bio.width.bi)
# For monofid piliform scale
rep.bio.width.mono <- rpt(formula = Largeur.Base ~ (1|SpeSpot), grname = c("SpeSpot"), data = ScalesMono, datatype = "Gaussian", nboot = 1000, nperm= 1000)
summary(rep.bio.width.mono)


### Calculation of median coefficient of variation for SpeSpot
## For scale density
density.bio.CV <- as.data.frame(matrix(NA, nrow = 6, ncol = 5 ))
names(density.bio.CV) <- c("SpeSpot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(density.bio.CV)){
  spespot <- levels(Density$SpeSpot)[i]
  MEAN <- mean(Density$densite.pedicelle[which(Density$SpeSpot == spespot)])
  SD <- sd(Density$densite.pedicelle[which(Density$SpeSpot == spespot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  density.bio.CV$speSpot[i] <- spespot
  density.bio.CV$mean[i] <- MEAN
  density.bio.CV$sd[i] <- SD
  density.bio.CV$CV[i] <- CV
  density.bio.CV$RSD[i] <- RSD*100
}
summary(density.bio.CV$RSD)

## For scale length
# For bifid piliform scale
length.biS.bio.CV <- as.data.frame(matrix(NA, nrow = 6, ncol = 5 ))
names(length.biS.bio.CV) <- c("SpeSpot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(length.biS.bio.CV)){
  spespot <- levels(ScalesBi$SpeSpot)[i]
  MEAN <- mean(ScalesBi$LongueurTot1[which(ScalesBi$SpeSpot == spespot)])
  SD <- sd(ScalesBi$LongueurTot1[which(ScalesBi$SpeSpot == spespot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  length.biS.bio.CV$speSpot[i] <- spespot
  length.biS.bio.CV$mean[i] <- MEAN
  length.biS.bio.CV$sd[i] <- SD
  length.biS.bio.CV$CV[i] <- CV
  length.biS.bio.CV$RSD[i] <- RSD*100
}
summary(length.biS.bio.CV$RSD)

# For monofid piliform scale
length.monoS.bio.CV <- as.data.frame(matrix(NA, nrow = 6, ncol = 5 ))
names(length.monoS.bio.CV) <- c("SpeSpot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(length.monoS.bio.CV)){
  spespot <- levels(ScalesMono$SpeSpot)[i]
  MEAN <- mean(ScalesMono$LongueurTot1[which(ScalesMono$SpeSpot == spespot)])
  SD <- sd(ScalesMono$LongueurTot1[which(ScalesMono$SpeSpot == spespot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  length.monoS.bio.CV$speSpot[i] <- spespot
  length.monoS.bio.CV$mean[i] <- MEAN
  length.monoS.bio.CV$sd[i] <- SD
  length.monoS.bio.CV$CV[i] <- CV
  length.monoS.bio.CV$RSD[i] <- RSD*100
}
summary(length.monoS.bio.CV$RSD)

## For scale wifth
# For bifid piliform scale
width.biS.bio.CV <- as.data.frame(matrix(NA, nrow = 6, ncol = 5 ))
names(width.biS.bio.CV) <- c("SpeSpot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(width.biS.bio.CV)){
  spespot <- levels(ScalesBi$SpeSpot)[i]
  MEAN <- mean(ScalesBi$Largeur.Base[which(ScalesBi$SpeSpot == spespot)])
  SD <- sd(ScalesBi$Largeur.Base[which(ScalesBi$SpeSpot == spespot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  width.biS.bio.CV$speSpot[i] <- spespot
  width.biS.bio.CV$mean[i] <- MEAN
  width.biS.bio.CV$sd[i] <- SD
  width.biS.bio.CV$CV[i] <- CV
  width.biS.bio.CV$RSD[i] <- RSD*100
}
summary(width.biS.bio.CV$RSD)

# For monofid piliform scale
width.monoS.bio.CV <- as.data.frame(matrix(NA, nrow = 6, ncol = 5 ))
names(width.monoS.bio.CV) <- c("SpeSpot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(width.monoS.bio.CV)){
  spespot <- levels(ScalesMono$SpeSpot)[i]
  MEAN <- mean(ScalesMono$Largeur.Base[which(ScalesMono$SpeSpot == spespot)])
  SD <- sd(ScalesMono$Largeur.Base[which(ScalesMono$SpeSpot == spespot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  width.monoS.bio.CV$speSpot[i] <- spespot
  width.monoS.bio.CV$mean[i] <- MEAN
  width.monoS.bio.CV$sd[i] <- SD
  width.monoS.bio.CV$CV[i] <- CV
  width.monoS.bio.CV$RSD[i] <- RSD*100
}
summary(width.monoS.bio.CV$RSD)


# --------------------------------------------------------
# Technical repeatability - measure repeatability
# --------------------------------------------------------
# We create a variable SpeSpot which is the interaction of species and spot on the forewing
Density$Spot <- interaction(Density$Code, Density$Zone)
Density$SpeSex <- interaction(Density$Genus, Density$sex)
ScalesBi$Spot <- interaction(ScalesBi$Code, ScalesBi$Veine)
ScalesMono$Spot <- interaction(ScalesMono$Code, ScalesMono$Veine)

### Calculation of the value of R with rpt function
# We look at the distribution of scale density, scale length and scale width
hist(Density$densite.pedicelle, nclass = 30)
hist(ScalesBi$LongueurTot1, nclass = 30)
hist(ScalesBi$Largeur.Base, nclass = 30)
hist(ScalesMono$LongueurTot1, nclass = 30)
hist(ScalesMono$Largeur.Base, nclass = 30)
# We consider that the variable are normally distributed because it is at least the closest distribution available in rptR
## For scale density
rep.tech.density <- rpt(formula = densite.pedicelle ~ Genus + SpeSex + (1|Code) + (1|Spot), grname = c("Code", "Spot"), data = Density, datatype = "Gaussian", nboot = 1000, nperm = 1000, adjusted = FALSE)
summary(rep.tech.density)
## For scale length
# For bifid piliform scale
rep.tech.length.bi <- rpt(formula = LongueurTot1 ~ Genus + sex + (1|Code) + (1|Spot), grname = c("Code", "Spot"), data = ScalesBi, datatype = "Gaussian", nboot = 1000, nperm= 1000, adjusted = FALSE)
summary(rep.tech.length.bi)
# For monofid piliform scale
rep.tech.length.mono <- rpt(formula = LongueurTot1 ~ Genus + sex + (1|Code) + (1|Spot), grname = c("Code", "Spot"), data = ScalesMono, datatype = "Gaussian", nboot = 1000, nperm= 1000, adjusted = FALSE)
summary(rep.tech.length.mono)
## For scale width
# For bifid piliform scale
rep.tech.width.bi <- rpt(formula = Largeur.Base ~ Genus + sex + (1|Code) + (1|Spot), grname = c("Code", "Spot"), data = ScalesBi, datatype = "Gaussian", nboot = 1000, nperm= 1000, adjusted = FALSE)
summary(rep.tech.width.bi)
# For monofid piliform scale
rep.tech.width.mono <- rpt(formula = Largeur.Base ~ Genus + sex + (1|Code) + (1|Spot), grname = c("Code", "Spot"), data = ScalesMono, datatype = "Gaussian", nboot = 1000, nperm= 1000, adjusted = FALSE)
summary(rep.tech.width.mono)


### Calculation of median coefficient of variation for SpeSpot
## For scale density
density.tech.CV <- as.data.frame(matrix(NA, nrow = 60, ncol = 5 ))
names(density.tech.CV) <- c("Spot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(density.tech.CV)){
  spot <- levels(Density$Spot)[i]
  MEAN <- mean(Density$densite.pedicelle[which(Density$Spot == spot)])
  SD <- sd(Density$densite.pedicelle[which(Density$Spot == spot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  density.tech.CV$Spot[i] <- spot
  density.tech.CV$mean[i] <- MEAN
  density.tech.CV$sd[i] <- SD
  density.tech.CV$CV[i] <- CV
  density.tech.CV$RSD[i] <- RSD*100
}
summary(density.tech.CV$RSD)

## For scale length
# For bifid piliform scale
length.biS.tech.CV <- as.data.frame(matrix(NA, nrow = 60, ncol = 5 ))
names(length.biS.tech.CV) <- c("Spot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(length.biS.tech.CV)){
  spot <- levels(ScalesBi$Spot)[i]
  MEAN <- mean(ScalesBi$LongueurTot1[which(ScalesBi$Spot == spot)])
  SD <- sd(ScalesBi$LongueurTot1[which(ScalesBi$Spot == spot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  length.biS.tech.CV$speSpot[i] <- spot
  length.biS.tech.CV$mean[i] <- MEAN
  length.biS.tech.CV$sd[i] <- SD
  length.biS.tech.CV$CV[i] <- CV
  length.biS.tech.CV$RSD[i] <- RSD*100
}
summary(length.biS.tech.CV$RSD)

# For monofid piliform scale
length.monoS.tech.CV <- as.data.frame(matrix(NA, nrow = 60, ncol = 5 ))
names(length.monoS.tech.CV) <- c("Spot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(length.monoS.tech.CV)){
  spot <- levels(ScalesMono$Spot)[i]
  MEAN <- mean(ScalesMono$LongueurTot1[which(ScalesMono$Spot == spot)])
  SD <- sd(ScalesMono$LongueurTot1[which(ScalesMono$Spot == spot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  length.monoS.tech.CV$speSpot[i] <- spot
  length.monoS.tech.CV$mean[i] <- MEAN
  length.monoS.tech.CV$sd[i] <- SD
  length.monoS.tech.CV$CV[i] <- CV
  length.monoS.tech.CV$RSD[i] <- RSD*100
}
summary(length.monoS.tech.CV$RSD)

## For scale wifth
# For bifid piliform scale
width.biS.tech.CV <- as.data.frame(matrix(NA, nrow = 60, ncol = 5 ))
names(width.biS.tech.CV) <- c("Spot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(width.biS.tech.CV)){
  spot <- levels(ScalesBi$Spot)[i]
  MEAN <- mean(ScalesBi$Largeur.Base[which(ScalesBi$Spot == spot)])
  SD <- sd(ScalesBi$Largeur.Base[which(ScalesBi$Spot == spot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  width.biS.tech.CV$speSpot[i] <- spot
  width.biS.tech.CV$mean[i] <- MEAN
  width.biS.tech.CV$sd[i] <- SD
  width.biS.tech.CV$CV[i] <- CV
  width.biS.tech.CV$RSD[i] <- RSD*100
}
summary(width.biS.tech.CV$RSD)

# For monofid piliform scale
width.monoS.tech.CV <- as.data.frame(matrix(NA, nrow = 60, ncol = 5 ))
names(width.monoS.tech.CV) <- c("Spot", "mean", "sd", "CV", "RSD")
for (i in 1:nrow(width.monoS.tech.CV)){
  spot <- levels(ScalesMono$Spot)[i]
  MEAN <- mean(ScalesMono$Largeur.Base[which(ScalesMono$Spot == spot)])
  SD <- sd(ScalesMono$Largeur.Base[which(ScalesMono$Spot == spot)])
  CV <- SD/MEAN
  RSD <- abs(CV)
  width.monoS.tech.CV$speSpot[i] <- spot
  width.monoS.tech.CV$mean[i] <- MEAN
  width.monoS.tech.CV$sd[i] <- SD
  width.monoS.tech.CV$CV[i] <- CV
  width.monoS.tech.CV$RSD[i] <- RSD*100
}
summary(width.monoS.tech.CV$RSD)



#_____________________________________________________________________________________________________
#### Supplementary result 1 - Link between physical and biologically relevant descriptors of transparency. ####
#_____________________________________________________________________________________________________

# -------------------------------------------------
# Data, libraries and home-made function needed 
# -------------------------------------------------
library(pavo)
library(nlme)
library(car)
library(caper)
getColSpace<-function(datatable, vision.mdl=c("VS", "UVS", "HUM"), illum=c("ws", "fs", "el", "lg"), conformation=c("transmission", "reflexion"), qcatch=c("Qi", "fi"), relative=c(FALSE, TRUE), colorspace){
  vs = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/sensitivity-curves/VS_shearwater.txt', header = T)
  colnames(vs) = c('wl','u','s','m','l','q')
  uvs = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/sensitivity-curves/UVS_bluetit.txt', header = T)
  colnames(uvs) = c('wl','u','s','m','l','q')
  illfs = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/light-environments/LuxFS.txt', header = T)
  ill_fs = illfs[,2]
  illel = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/light-environments/LuxEL.txt', header = T)
  ill_el = illel[,2]
  illws = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/light-environments/LuxWS.txt', header = T)
  ill_ws = illws[,2]
  illlg = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/light-environments/LuxLG.txt', header = T)
  ill_lg = illlg[,2]
  fondfeuil = read.delim('C:/Users/Charline/Documents/These_Transparency-mimicry/Data_raw/Vision_modeling/backgrounds/FondCGfeuilles.txt', header = T)
  fond_feuille = fondfeuil[,2]
  
  if (vision.mdl=="UVS") {
    vis<-uvs 
    nvis<-c(1, 1.9, 2.7, 2.7)
    webvis<-0.1
    vision<-"uvs" #pour LWS
    visdf = vis[,c(1:5)]
    achdf = vis[,6]
    webachr=0.2
  }
  if (vision.mdl=="VS") {
    vis<-vs
    nvis<-c(1, 0.7, 1, 1.4)
    vision<-"vs"
    webvis<-0.1 #pour LWS
    visdf = vis[,c(1:5)]
    achdf = vis[,6]
    webachr=0.2
  }
  if (vision.mdl=="HUM")  {
    vis<-hum
    nvis<-c(1,16,32)
    vision<-"hum"
    webvis<-0.018 #pour LWS
    visdf = vis[,c(1:4)]
    achdf = vis[,5]
    webachr=0.11 
  }
  
  if (conformation=="transmission") {
    data<-datatable
  }
  if (conformation=="reflexion") {
    data<-as.data.frame(matrix(NA, nrow=dim(datatable)[1], ncol=dim(datatable)[2]))
    names(data)<-names(datatable)
    for (j in 2:dim(datatable)[2]){
      for (i in 1:dim(datatable)[1]){
        data[i,j]<-datatable[i,j]*datatable[i,j]*fond_feuille[i]/10000
      }
    }
    data[,1]<-datatable$wl
  }
  
  if (illum=="ws"){ 
    illdf<-ill_ws
  }
  if (illum=="fs"){ 
    illdf<-ill_fs
  }
  if (illum=="el"){ 
    illdf<-ill_el
  }
  if (illum=="lg"){ 
    illdf<-ill_lg
  }
  
  data.vismodel<-vismodel(data, visual = visdf, achromatic = achdf, illum = illdf, trans=1, qcatch=qcatch,relative=relative)
  data.colspace<-colspace(data.vismodel, space = colorspace, qcatch=qcatch)
  
  return(data.colspace)  
}
# phylogenetic tree of the species from the study
tree <- read.nexus("C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/Documents_pour-article/arbreM2_mai2020_simple.trees")
# Dataframe with all spectra
specall <- read.delim(file = "C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/Lepidoptera_Transparenct-Mimicry_SpecAll.txt", sep = "\t", dec = ".", check.names = FALSE)
specall <- as.rspec(specall)
# Information about spectra and specimens
specinfos <- read.delim(file = "C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/Lepidoptera_Transparenct-Mimicry_SpecInfos.txt", header = TRUE, sep ="\t", dec = ".")
# Exploitation of raw data
# We extract mean transmittance over 300-700nm
specall.shape <- summary(specall, subset = "B2")
specall.shape$SpectrumID <- row.names(specall.shape)
# We extract the xyz coordinates in tetrahedral space and luminance from a visual model with UVs vision and large gap ambiant light
specall.colspace <- getColSpace(datatable = specall, vision.mdl = "UVS", illum = "lg", conformation = "transmission", qcatch = "Qi", relative = TRUE, colorspace = "tcs")
specall.colspace$SpectrumID <- row.names(specall.colspace)
# We merge the different dataframe to get the dataset to test for the link between mean transmittance over 300-700nm (B2) and biologically relevent descriptors (xyzL)
linkB2xyzL.data <- merge(merge(specall.shape, specall.colspace[,c("SpectrumID", "x", "y", "z", "lum")], by = "SpectrumID"), specinfos, by = "SpectrumID")

# -------------------------------------------------
# Linear mixed model 
# -------------------------------------------------
linkB2xyzL.data[,c("B2", "x", "y", "z", "lum")] <- scale(linkB2xyzL.data[,c("B2", "x", "y", "z", "lum")])
linkB2xyzL.mdl.lme<-lme(fixed = B2 ~ x + y + z + lum, data = linkB2xyzL.data, random = ~ 1|SpecimenID, method = "ML")
#results of the linear mixed model (table 1a)
summary(linkB2xyzL.mdl.lme)
#Table of deviance (table 1b)
Anova(linkB2xyzL.mdl.lme)


# -------------------------------------------------
# PGLS model
# -------------------------------------------------
# We create a dataframe with mean value of x, y, z, lum and B2 for each species
linkB2xyzL.data.mean <- as.data.frame(matrix(NA, nrow = length(unique(linkB2xyzL.data$SpecimenID)), ncol = 11))
names(linkB2xyzL.data.mean) <- c("SpecimenID", "SpeShort", "Tip.Label", "Genus", "species", "mimicry.ring", "B2", "x", "y", "z", "lum")
for (i in 1:length(unique(linkB2xyzL.data$SpecimenID))){
  specimen <- unique(linkB2xyzL.data$SpecimenID)[i]
  mimicry <- unique(linkB2xyzL.data$mimicry.ring[which(linkB2xyzL.data$SpecimenID == specimen)])
  if (mimicry != "agnosia") toMean <- which(linkB2xyzL.data$SpecimenID == specimen) else toMean <- which(linkB2xyzL.data$SpecimenID == specimen & linkB2xyzL.data$Spot != "s1") 
  linkB2xyzL.data.mean$SpecimenID[i] <- specimen
  linkB2xyzL.data.mean$SpeShort[i] <- unique(linkB2xyzL.data$SpeShort[which(linkB2xyzL.data$SpecimenID == specimen)])
  linkB2xyzL.data.mean$Tip.Label[i] <- unique(linkB2xyzL.data$Tip.Label[which(linkB2xyzL.data$SpecimenID == specimen)])
  linkB2xyzL.data.mean$Genus[i] <- unique(linkB2xyzL.data$Genus[which(linkB2xyzL.data$SpecimenID == specimen)])
  linkB2xyzL.data.mean$species[i] <- unique(linkB2xyzL.data$species[which(linkB2xyzL.data$SpecimenID == specimen)])
  linkB2xyzL.data.mean$mimicry.ring[i] <- mimicry
  linkB2xyzL.data.mean$B2[i] <- mean(linkB2xyzL.data$B2[toMean])
  linkB2xyzL.data.mean$x[i] <- mean(linkB2xyzL.data$x[toMean])
  linkB2xyzL.data.mean$y[i] <- mean(linkB2xyzL.data$y[toMean])
  linkB2xyzL.data.mean$z[i] <- mean(linkB2xyzL.data$z[toMean])
  linkB2xyzL.data.mean$lum[i] <- mean(linkB2xyzL.data$lum[toMean])
}

# We create a comparative.data object to run pgls model with caper
linkB2xyzL.compdata <- comparative.data(phy = tree, data = linkB2xyzL.data.mean, names.col = Tip.Label, vcv = TRUE)
linkB2xyzL.compdata$data[,c("B2", "x", "y", "z", "lum")]<-scale(linkB2xyzL.compdata$data[,c("B2", "x", "y", "z", "lum")])
# We compute pgls model
linkB2xyzL.pgls <- pgls(formula = B2 ~ x + y + z + lum, data=linkB2xyzL.compdata, lambda = "ML")
# Results of the pgls model(table 1c)
summary(linkB2xyzL.pgls)


#_____________________________________________________________________________________________________
#### Supplementary result 2 - Link between nanostructure density and nanostructure type ####
#_____________________________________________________________________________________________________

# -------------------------------------------------
# Data, libraries and home-made function needed 
# -------------------------------------------------
library(phytools)
# Datatable with all the structural features and mean transmittance
DataTot <- read.delim(file = "C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/article_Mimicry-M2Charline/Data/Lepidoptera_Transparenct-Mimicry_Structure-Transmission-Infos-all.txt", header = TRUE, sep = "\t", dec = ".")
# phylogenetic tree of the species from the study
tree <- read.nexus("C:/Users/Charline/Documents/These_Transparency-mimicry/Article/articleM2_mimetisme-transparence/Documents_pour-article/arbreM2_mai2020_simple.trees")

# -------------------------------------------------
# Analyses - phylogenetic anova
# -------------------------------------------------
NanoType <- DataTot$NanoType
names(NanoType) <- DataTot$TipLabel
NanoDens <- DataTot$NanoDensity.mean
names(NanoDens) <- DataTot$TipLabel
nano.phyaov <- phylANOVA(tree = tree, x = NanoType, y = NanoDens, nsim = 1000, posthoc = TRUE, p.adj = "bonferroni")
# Results from the phylogenetic anova (table 2a and 2b)
nano.phyaov
