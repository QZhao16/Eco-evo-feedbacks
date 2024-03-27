# Set path
Root_Path <- 'C:/zhao088/Desk20210304/Manuscript.eco-evo/6 Lake Oneida/MDR S-map'
OPath <- paste(Root_Path,sep='')
setwd(OPath)


#First <- function(){ .libPaths("C:/Program Files/R/R-3.5.0/library") }
#First()
#gc()
#rm()




#####################################################
library('rEDM')
library('dplyr')
library('tidyr')
library('igraph')
library('MASS')
library('NetIndices')
library('Cairo')
library('vegan')
library('Kendall') # Kendall's tau test for the convergence of CCM

library('doParallel')
library('parallel')
library('foreach')
library('glmnet')
library('psych') #trace
library('pracma')


# unit (/Littre)
  new_alg_zo<-read.csv("i.phyto.zoo.abu.biom.seas.csv")
  new_alg_zo<-as.data.frame(new_alg_zo)# having year
  dim(new_alg_zo)
  # scale
  alga_zo<-apply(new_alg_zo[,2:dim(new_alg_zo)[[2]]], 2, scale) # no year
  dim(alga_zo) # 161  27
  head(alga_zo) # 161 19
  colnames(alga_zo)
  
  E_alg_zoo<-NULL
  # change 1
  for (i in 1:dim(alga_zo)[[2]] ) {
    
    sim_r <- simplex(alga_zo[,i],lib=c(1,floor(nrow(alga_zo)/2)),pred=c(floor(nrow(alga_zo)/2)+1,nrow(alga_zo)),E=c(2:9))
    #The optimal embedding dimension determined by maximizing rho 
    (E_r <-sim_r[which.min(sim_r$mae),"E"][1]) 
    E_alg_zoo <- rbind(E_alg_zoo,E_r)}
  
  dim(E_alg_zoo)  
  #
  write.csv(E_alg_zoo,"E_alg_zoo(mae24).csv")
  
  
  ##############################################
  # part 6          causality+time delay       #
  ##############################################
  #
  E_alg_zoo<-read.csv("E_alg_zoo(mae24).csv")
  E_alg_zoo[,2]
  length(E_alg_zoo[,2]) # 23
  
  ########## ccm: 19 speceis
  i_j_rho<-NULL
  j_rho<-NULL
  
  for (i in 1:ncol(alga_zo) ) {
    #195=dim(alga_zo)[[2]]
    for (j in 1:ncol(alga_zo) ) {
      #remov(i)
      lag_rho<-NULL
      
      for (lags in -2:0) {   
        E_star<- E_alg_zoo[i,2]
        
        dat_temp <- alga_zo[,c(i,j)]  # in ccm, use i to estimate j, so use E of i; but in practice we do use j to estimate i
        
        block_temp = data.frame(laf(dat_temp[,1],dat_temp[,2],lagf=lags))  ## size=dat_temp[,2]
        colnames(block_temp) = c("RUE","Diversity")                        ## laf() function here is to just change clolumn names
        
        ##lib_ccm <- c(1,NROW(block_temp))  old cold
        lib_ccm <- c(1,NROW(block_temp))
        ##pred_ccm <- make_pred_nozero(alga_zo[,i],E_star)   old cold
        pred_ccm <- make_pred_nozero(block_temp$RUE,E_star)
        ##
        df.out.ccm <- ccm(block=block_temp,
                          E=E_star,
                          lib=lib_ccm,        ###[1]   1 329 (from 1 to end to attractor reconstruction)
                          pred = pred_ccm,    ###[1]   1 328
                          lib_column = 1,     ## speceis 1
                          target_column = 2,  ## speceis i
                          lib_sizes = NROW(block_temp),    #lib_sizes = NROW(block_temp),
                          exclusion_radius=0,
                          random_libs = FALSE,
                          num_sample=1,
                          tp = 0)
        
        lag_rho<-cbind(lag_rho, df.out.ccm$rho)} 
      
      j_rho<-cbind(i,j,lag_rho)
      i_j_rho<-rbind(i_j_rho,j_rho)
    }
    
  }
  
  # do not delete
  write.csv(i_j_rho,"rho_realdata.csv") # best lag+ largest rho
  
  
  ####################################################################
  ######################### find maximum lag #########################
  ####################################################################
  
  rho_realdata<-read.csv("rho_realdata.csv")
  head(rho_realdata)
  dim(rho_realdata) #256 16
  rho_realdata[,4]
  
  ##
  a<-NULL
  b<-NULL
  c<-NULL
  d<-NULL
  
  out <- data.frame(i=rep(NA,dim(rho_realdata)[[1]]),j=rep(NA,dim(rho_realdata)[[1]]),lag=rep(NA,dim(rho_realdata)[[1]]), maximum_rho=rep(NA,dim(rho_realdata)[[1]]))
  
  for (i in 1:dim(rho_realdata)[[1]]) { 
    a<-rho_realdata[i,4:dim(rho_realdata)[[2]]]
    a[is.na(a)] <- 0
    a[which(a<0)]<-0
    b<-which.max(a)
    c<-a[which.max(a)]
    
    out[i,1]<- rho_realdata$i[i]
    out[i,2]<-rho_realdata$j[i]
    
    out[i,3]<- (b-3)
    out[i,4]<-c[[1]]
    
  }
  
  dim(out) #256   4
  head(out)
  
  #
  write.csv(out,"rho_realdata_pluslag.csv")
  
  
  ##################################################
  # unit (/Littre)
  new_alg_zo<-read.csv("i.phyto.zoo.abu.biom.seas.csv")
  new_alg_zo<-as.data.frame(new_alg_zo)# having year
  dim(new_alg_zo)
  # scale
  alga_zo<-apply(new_alg_zo[,2:dim(new_alg_zo)[[2]]], 2, scale) # no year
  dim(alga_zo) # 161  27
  head(alga_zo) # 161 19
  colnames(alga_zo)
  
  #
  maximum_lag_realdata<-read.csv("rho_realdata_pluslag.csv")
  head(maximum_lag_realdata)
  length(maximum_lag_realdata[,2])  #361
  maximum_lag_realdata[,4]          #361
  
  #
  E_alg_zoo<-read.csv("E_alg_zoo(mae24).csv")
  E_alg_zoo[,2]
  
  ########################## rho_difference between maxmum-rho and minmum-rho
  i_j_rho<-NULL
  j_rho<-NULL
  
  #
  for (i in 1:dim(alga_zo)[[2]]) {
    #16=dim(alga_zo)[[2]]
    for (j in 1:dim(alga_zo)[[2]]) {
      
      rho_difference<-NULL
      
      E_star<- E_alg_zoo[i,2]
      k<-(i-1)*dim(alga_zo)[[2]]+j
      
      dat_temp <- alga_zo[,c(i,j)]
      
      block_temp = data.frame(laf(dat_temp[,1],dat_temp[,2],lagf=maximum_lag_realdata[k,4])) # fro -lag(time daly) 
      colnames(block_temp) = c("RUE","Diversity")
      
      out.temp <- ccm(block=block_temp,
                      E=E_star,
                      lib_column = 1,
                      target_column = 2,
                      lib_sizes = seq(9, nrow(block_temp), by = 10),
                      random_libs = F, 
                      tp=0) # no lag,because block_temp has already lagged
      #num_samples = 100, 
      #random_libs = F, 
      #replace = TRUE)
      
      #out.temp$rho[which(out.temp$rho<0)]<-0              
      out.temp$rho[is.na(out.temp$rho)]<-0                
      out.temp$rho[which(!is.finite(out.temp$rho))] <- 0  
      
      out.temp_means <- ccm_means(out.temp)
      differenc<-out.temp_means$rho[length(out.temp_means$rho)]-out.temp_means$rho[1]
      rho_difference<-cbind(rho_difference, differenc) 
      
      # Kendall's ?? test and Student's t-test aresiginificant
      np<-nrow(block_temp);
      crirho <- qt(0.95,np-1)/(np-2+qt(0.95,np-1)^2);
      kend <- MannKendall(out.temp$rho);  # Kendall's tau test for mononotic increase,both (1) Kendall's tau and (2) terminal rho are significantly larger than zero
      significance_mononotic <- (kend$tau[1]>0)*(kend$sl[1]<0.05)*(out.temp_means$rho[length(out.temp_means$rho)]>crirho) # ccm.sig records the significance of each CCM
      P_value_mono<-kend$sl[1]
      rho_mono<-out.temp_means$rho[length(out.temp_means$rho)]
      
      j_rho<-cbind(i,j,maximum_lag_realdata[k,4],rho_difference,significance_mononotic,P_value_mono,rho_mono)
      #j_rho<-cbind(i,j,maximum_lag_realdata[k,4],rho_difference)
      i_j_rho<-rbind(i_j_rho,j_rho)
      
    }
  }
  
  write.csv(i_j_rho,"rho_difference_E.csv")
  
  
  
  
  ## cbind  maximum_lag + rho_difference
  maximum_lag_realdata<-read.csv("rho_realdata_pluslag.csv")
  length(maximum_lag_realdata[,2])  
  maximum_lag_realdata[,3]   
  # change 2.2
  rho_diff<-read.csv("rho_difference_E.csv")
  rho_diff$significance<-ifelse(rho_diff$significance_mononotic==1,1,0) #>0 as convergence
  #rho_diff$significance<-ifelse(rho_diff$differenc>0.0,1,0) #>0 as convergence
  head(rho_diff)
  
  maximum_lag_rho_difference<-cbind(rho_diff,maximum_lag_realdata[,5])
  write.csv(maximum_lag_rho_difference,"maxlag_rhodifer.csv")
  
  
  
  ###################################
  ################ceating surrogates
  # unit (/Littre)
  new_alg_zo<-read.csv("i.phyto.zoo.abu.biom.seas.csv")
  new_alg_zo<-as.data.frame(new_alg_zo)# having year
  dim(new_alg_zo)
  # scale
  alga_zo<-apply(new_alg_zo[,2:dim(new_alg_zo)[[2]]], 2, scale) # no year
  dim(alga_zo) # 161  27
  head(alga_zo) # 161 19
  colnames(alga_zo)
  
  
  #add date
  #Date<-seq(ISOdate(2000,1,4), by = "8 DSTdays", length.out = 695)
  #date <-read.csv("date.csv")
  new_alg<-read.csv("i.phyto.zoo.abu.biom.seas.csv")
  Date<-as.Date(new_alg$sampledate,"%m/%d/%Y")
  #Date<-Date[1:308]  # no 2017-year & start (08-May-1991)
  #as.Date(date$year_date, origin = "08-Jan-1992") 
  #write.csv(Date,"Date.csv")
  # load E
  E_alg_zoo<-read.csv("E_alg_zoo(mae24).csv")
  E_alg_zoo[,2]
  
  ######################## generating surrogates
  
  surrogate_allspeceis<-NULL
  for (j in 1:dim(alga_zo)[[2]]) {
    #remov(i)
    
    df.in<-data.frame (date = Date, Diversity = alga_zo[,j])
    #to surrogates
    
    set.seed(599213) 
    
    out_diversity <- yearday_anom(df.in$date,df.in$Diversity)
    Diversity.bar <- out_diversity$mean
    Diversity.tilde <- out_diversity$anomaly
    
    
    # Diversity surrogates  
    Diversity.surrs <- do.call(cbind,
                               lapply(1:100, function(ii) {
                                 I_na <- is.na(Diversity.tilde)
                                 out <- Diversity.bar
                                 out[I_na] <- NA
                                 out[!I_na] <- out[!I_na] + sample(Diversity.tilde[!I_na],sum(!I_na),replace = FALSE)
                                 return(out)
                               }))
    
    surrogate_allspeceis<-cbind(surrogate_allspeceis,Diversity.surrs)
  }
  
  dim(surrogate_allspeceis) #303 1500
  # do not delte
  #write.csv(surrogate_allspeceis,"surrogate_allspeceis.csv")
  save(surrogate_allspeceis, file = paste('surrogate_allspeceis','.RData',sep = ""))
  
  load("surrogate_allspeceis.RData")
  dim(surrogate_allspeceis) #300 1900
  # try to overlap ture+surrogate
  plot(surrogate_allspeceis[,1], type="l")
  lines(alga_zo[,1],col="red")
  
  #################################################################################
  load("surrogate_allspeceis.RData")
  dim(surrogate_allspeceis) #300 1900
  
  order_lag_causality<-read.csv("maxlag_rhodifer.csv")
  order_lag_causality<-subset(order_lag_causality,significance==1)  # important:select rho_difference>0.1
  dim(order_lag_causality) # 251   8
  
  E_alg_zoo<-read.csv("E_alg_zoo(mae24).csv")
  E_alg_zoo[,2]
  
  
  ###################### rho of surrogates: Using real VS surrogates ##
  
  i_j_rho<-NULL
  j_rho<-NULL
  
  for (i in 1:length(order_lag_causality$i)) {
    
    lag_rho<-NULL
    
    E_star<- E_alg_zoo[order_lag_causality$i[i],2]
    max_range<-c(1:100)+(order_lag_causality$j[i]-1)*100
    #k<-(i-1)*195+j
    
    for (i_surr in  max_range){
      dat_temp <- data.frame(RUE=alga_zo[,order_lag_causality$i[i]], Diversity=surrogate_allspeceis[,i_surr])
      
      lag_rho_surroga <- do.call(cbind,
                                 lapply(0, function(ii) {
                                   block_temp = data.frame(laf(dat_temp[,1],dat_temp[,2],lagf=ii))  
                                   colnames(block_temp) = c("RUE","Diversity")                        
                                   lib_ccm <- c(1,NROW(block_temp))
                                   pred_ccm <- make_pred_nozero(block_temp$RUE,E_star)
                                   
                                   df.out.ccm <- ccm(block=block_temp,
                                                     E=E_star,
                                                     lib=lib_ccm,        ###[1]   1 329 (from 1 to end to attractor reconstruction)
                                                     pred = pred_ccm,    ###[1]   1 328
                                                     lib_column = 1,     ## speceis 1
                                                     target_column = 2,  ## speceis i
                                                     lib_sizes = NROW(block_temp),    #lib_sizes = NROW(block_temp),
                                                     exclusion_radius=0,
                                                     random_libs = FALSE,
                                                     num_sample=1,
                                                     tp = 0)
                                   return(df.out.ccm$rho)
                                 }))
      
      out.temp<-lag_rho_surroga[which.max(lag_rho_surroga)]
      lag_rho<-cbind(lag_rho, out.temp)} 
    
    j_rho<-cbind(order_lag_causality$i[i],order_lag_causality$j[i],lag_rho)
    i_j_rho<-rbind(i_j_rho,j_rho)
    
    
  }
  
  save(i_j_rho, file = paste('rho_in_surrodata(mae24)','.RData',sep = ""))
  write.csv(i_j_rho, "rho_in_surrodata(mae24).csv")
  
  load("rho_in_surrodata(mae24).RData")
  dim(i_j_rho) #232 102
  
  

  # change 3
  ######
  order_lag_causality<-read.csv("maxlag_rhodifer.csv")
  order_lag_causality<-subset(order_lag_causality,significance==1)  #select rho_difference>0.0, then to 
  dim(order_lag_causality) # 
  order_lag_causality$maximum_lag_realdata
  colnames(order_lag_causality)<- c("number","X",	"i","j", "lag",	"differenc","significance_mononotic","P_value_mono","rho_mono", "significance",	"maximum_lag_realdata")
  
  #
  rho_surrodata1<-read.csv("rho_in_surrodata(mae24).csv") #these surrogates only from rho_difference>0.0
  rho_surrodata1[,4]
  rho_surrodata<-rho_surrodata1[1:dim(rho_surrodata1)[[1]],4:dim(rho_surrodata1)[[2]]]
  dim(rho_surrodata) # 
  colnames(rho_surrodata1)<-c("number", "i",	"j",	rep("out.temp", times=100) )
  
  ############# merge significance: real data (rho difference>0.1)+ surrogate data (max rho in real data> rho in surrogates)
  aa<-NULL
  bb<-NULL
  
  length(order_lag_causality$maximum_lag_realdata)
  # change 4
  ou <- data.frame(i=rho_surrodata1$i, j=rho_surrodata1$j, lag=order_lag_causality$lag, rho_difference=order_lag_causality$differenc, 
                   significance_mononotic=order_lag_causality$significance_mononotic,P_value_mono=order_lag_causality$P_value_mono,
                   rho_mono=order_lag_causality$rho_mono,significance=order_lag_causality$significance,
                   max.rho_fromlag_realdata=order_lag_causality$maximum_lag_realdata, significance_surrogates=rep(NA,dim(rho_surrodata1)[[1]]))
  
  for (ii in 1:dim(rho_surrodata)[[1]]) { 
    aa<-rho_surrodata[ii,]
    aa[is.na(aa)] <- 0
    aa[which(aa<0)]<-0
    #calculating p-values
    bb<-(sum(ou$max.rho_fromlag_realdata[ii] < aa) + 1) / (length(aa) + 1)
    ou$significance_surrogates[ii]<-bb 
  }
  
  dim(ou) #
  head(ou)
  
  write.csv(ou,"interaction_significance.csv") 
  
  
  
  ##################################
  ### In s-map, we do not need E ###
  ##################################
  library('rEDM')
  library('dplyr')
  library('tidyr')
  
  # unit (/Littre)
  new_alg_zo<-read.csv("i.phyto.zoo.abu.biom.seas.csv")
  new_alg_zo<-as.data.frame(new_alg_zo)# having year
  dim(new_alg_zo)
  # scale
  alga_zo<-apply(new_alg_zo[,2:dim(new_alg_zo)[[2]]], 2, scale) # no year
  dim(alga_zo) # 161  27
  head(alga_zo) # 161 19
  colnames(alga_zo)
  
  
  ### only extract causal interactions(rho>0) + order  
  
  interaction_significance<-read.csv("interaction_significance.csv")
  # change 5
  #out<-subset(interaction_significance,i!=j & significance_surrogates<0.05) 
  out<-subset(interaction_significance, significance_surrogates<0.05) 
  head(out)               
  dim(out) #670   7
  
  ###order it                
  order_lag_causality<-NULL
  for (ii in unique(out$i)) {
    # extract causality for each ii
    sub_lag_causality<-subset(out, i==ii) #>0.1
    # order ascending
    sub_lag_causality<-sub_lag_causality[order(sub_lag_causality$lag),] # ascending with lag
    # row bind
    order_lag_causality<- rbind(order_lag_causality,sub_lag_causality)
  }
  dim(order_lag_causality) #187   7
  
  
  # then merge E+links
  #############################################
  E_alg_zoo<-read.csv("E_alg_zoo(mae24).csv")
  E_alg_zoo[,2]
  
  for (ii in 1:length(order_lag_causality$i)) {order_lag_causality$Best_E[ii]<-E_alg_zoo[order_lag_causality$i[ii],2]
  order_lag_causality$total_links[ii]<-length(order_lag_causality$i[order_lag_causality$i==order_lag_causality$i[ii]])
  }
  dim(order_lag_causality) #[1] 114  9
  
  # then, order them
  ###################
  final_order_laga<-NULL
  
  for (ii in unique(order_lag_causality$i)) {
    # extract causality for each ii
    sub_order_lag_causality<-subset(order_lag_causality, i==ii)
    # order ascending
    sub_lag<-sub_order_lag_causality[order(sub_order_lag_causality$lag),] #ranking from negative lag to positive lag 
    # row bind
    final_order_laga<- rbind(final_order_laga,sub_lag)
  }
  
  write.csv(final_order_laga,"caus_ord_lag(mae).csv")
  dim(final_order_laga) #[1] 52  9
  
  
  #then
  #### significance: real+surrogates
  ####################### check mean+sd causal links/ links > E
  E_alg_zoo<-read.csv("E_alg_zoo(mae24).csv")
  E_alg_zoo[,2]
  
  ave<-read.csv("caus_ord_lag(mae).csv")
  dim(ave) #[1] 95 10
  length(unique(ave$i)) #[1] 15 species
  #count links > E
  length(which(ave$total_links>ave$Best_E)) #[1] 5
  
  ave[which(ave$total_links>ave$Best_E),"total_links"] #return total_links
  ave[which(ave$total_links>ave$Best_E),"Best_E"] #return Best_E
  ave[which(ave$total_links>ave$Best_E),"lag"]    #return lag
  
  sum(ave[which(ave$total_links>ave$Best_E),"total_links"]-ave[which(ave$total_links>ave$Best_E),"Best_E"])/length(ave[which(ave$total_links>ave$Best_E),"Best_E"]) 
  
  ################################
  
  

  
  ##------------------------ --------##
  ##----delete indirect links--------##
  ##------------------------ --------##
  library(igraph)
  library(MASS)
  library(NetIndices)
  require(Cairo)
  
  
  otago.links.data<-read.csv("caus_ord_lag(mae).csv")
  otago.links.data<-subset(otago.links.data, i!=j) # exclude i=j
  otago.links.data<-subset(otago.links.data, significance_surrogates<0.05) # p<0.05
  
  colnames(otago.links.data)
  otago.graph.c<-graph.edgelist(as.matrix(otago.links.data[,c(3,4)]))
  otago.adjmatrix.c<-get.adjacency(otago.graph.c,sparse=F)
  
  #
  diameter(otago.graph.c, unconnected=TRUE)
  #
  otago.adjmatrix.c[1,]
  otago.adjmatrix.c[2,]
  dim(otago.graph.c)
  otago.graph.c[]
  
  # otago.graph.c<-graph.edgelist(as.matrix(otago.links.data[,c(4,3)]))
  #
  
  
  ae=ae2=ae3=ae4=ae5=ae6=ae7=ae8=NULL
  ag=ah=ai=aj=ak=al=am=NULL
  for (i1 in 1:dim(otago.adjmatrix.c)[[2]] ) { 
    
    # via 1 transfer
    a<-otago.adjmatrix.c[i1,] #  1<-- 
    
    if ( length(which(a==1)) >=1 ) { 
      for (i2 in 1:length(which(a==1)) ) {  aa<- which(a==1)
      targe<-aa[i2]
      ab<-aa[-i2]
      for (i3 in ab ) { b<-otago.adjmatrix.c[, targe ] 
      ac<-which(b==1)
      if (i3 %in% ac ){ad<-cbind(targe,i3,i1)
      ae<- rbind(ae,ad)}    } 
      }
    }
    
    # via 2 transfer
    a<-otago.adjmatrix.c[i1,] #  1<-- 
    if ( length(which(a==1)) >=1 ) { 
      for (i2 in 1:length(which(a==1)) ) {  aa<- which(a==1)
      targe<-aa[i2]
      ab<-aa[-i2]
      for (i3 in ab ) { b<-otago.adjmatrix.c[, targe ] # 1th transfer
      ac<-which(b==1)
      ac<- ac[-which(ac==i1)] # remove i1, otherwise have closed loop
      
      treansfer2 <- do.call(rbind,
                            lapply( ac, function(targe2) {
                              c<-otago.adjmatrix.c[, targe2 ] # 2th transfer
                              af<-which(c==1)
                              if (i3 %in% af ){ag<-cbind(targe,targe2,i3,i1)
                              #ae2<- rbind(ae2,ag)
                              }
                              return(ag)
                            }))
      
      ae2<-rbind(ae2,treansfer2)    } 
      }
    }
    
    
    
  }
  
  
  # have vaule
  ae
  ae2

  
  
  
  ae<- as.data.frame(ae)
  if (ncol(ae) ==1) {ae<-t(ae) 
  #
  newmatrix<- as.data.frame( matrix(rep(0,times=ncol(ae)),nrow=1) )
  colnames(newmatrix)<-colnames(ae)
  for (i in 1:ncol(ae) ) { newmatrix[1,i]<- ae[1,i] }
  #
  ae<-c()
  ae<-newmatrix
  }
  
  
  ae2<- as.data.frame(ae2)
  if (ncol(ae2) ==1) {ae2<-t(ae2) 
  #
  newmatrix<- as.data.frame( matrix(rep(0,times=ncol(ae2)),nrow=1) )
  colnames(newmatrix)<-colnames(ae2)
  for (i in 1:ncol(ae2) ) { newmatrix[1,i]<- ae2[1,i] }
  #
  ae2<-c()
  ae2<-newmatrix
  }
  

  #============= first remove repeate ones within
  combi<- NULL
  if (nrow(ae)>0){
    for (i in dim(ae)[1] ) {  subse<- ae[i,]
    if ( length(subse) > length(unique(subse)) ) {combi<-cbind(combi,i)}} #ae<-ae[-i,]
    
    if (length(combi)>0){ae<-ae[-combi,]}
  }
  #
  combi<- NULL
  if (nrow(ae2)>0){
    for (i in dim(ae2)[1] ) {  subse<- ae2[i,]
    if ( length(subse) >  length(unique(subse)) ) {combi<-cbind(combi,i)} } #ae2<-ae2[-i,]
    
    if (length(combi)>0){ae2<-ae2[-combi,]}
  }
  #
  #============= second, remove repeated short path from long path
  ae
  ae2
  
  
  #=====
  
  # to the middle. I do not special care, only care (2)
  # (1) 3->8->1->13  & (2) 3->11->1
  
  # from the middle. I do not special care. just assort together as below
  # 3->11->1 & 19->3->17->1  
  
  # filling the middle. I do not special care. just assort together as below
  # 3->11->1 & 3->12->11->1
  
  
  ae<- as.data.frame(ae)
  if (ncol(ae) ==1) {ae<-t(ae) 
  #
  newmatrix<- as.data.frame( matrix(rep(0,times=ncol(ae)),nrow=1) )
  colnames(newmatrix)<-colnames(ae)
  for (i in 1:ncol(ae) ) { newmatrix[1,i]<- ae[1,i] }
  #
  ae<-c()
  ae<-newmatrix
  }
  
  
  ae2<- as.data.frame(ae2)
  if (ncol(ae2) ==1) {ae2<-t(ae2) 
  #
  newmatrix<- as.data.frame( matrix(rep(0,times=ncol(ae2)),nrow=1) )
  colnames(newmatrix)<-colnames(ae2)
  for (i in 1:ncol(ae2) ) { newmatrix[1,i]<- ae2[1,i] }
  #
  ae2<-c()
  ae2<-newmatrix
  }
  
  
  
  
  # change 5.2
  number_species=46
  need_delete=need_delet=need_deleted=NULL
  caus_ord<-read.csv("caus_ord_lag(mae).csv")
  
  
  for ( l in 1:number_species ) {   #l is end point 
    for ( k in 1:number_species ) { #k is start point
      
      #===== ae
      count_ae<-0
      
      if ( nrow(ae) == 0 ){
        count_ae<-0
        su_ae<-NULL
        su_ae<-as.data.frame(su_ae)
        nrow(su_ae)}
      
      if ( nrow(ae) != 0 ){
        
        sub_ae<- subset(ae,i1==l) 
        #
        if (nrow(sub_ae) ==0) {
          su_ae<-NULL
          su_ae<-as.data.frame(su_ae)
          nrow(su_ae)}
        #
        if (nrow(sub_ae) >=1) {
          su_ae<- subset(sub_ae,targe==k) 
          if ( nrow(su_ae) >=1 ) { sub_caus_ord<- subset(caus_ord, i==l) 
          su_caus_ord<- subset(sub_caus_ord, j==k) #compare with this (longest link)
          for(m in 1:nrow(su_ae) ) {
            
            compare_caus_ord<- subset(caus_ord,j==su_ae$targe[m])  # targe is start point
            compare_caus_or<- subset(compare_caus_ord,i==su_ae$i3[m])# i3 is end point, using this for next
            
            # compare with this
            if (compare_caus_or$max.rho_fromlag_realdata > su_caus_ord$max.rho_fromlag_realdata & compare_caus_or$lag > su_caus_ord$lag) {count_ae<- count_ae+1} } }
        }}
      
      
      #===== ae2
      count_ae2<-0
      
      if ( nrow(ae2) == 0 ){
        count_ae2<-0
        su_ae2<-NULL
        su_ae2<-as.data.frame(su_ae2)
        nrow(su_ae2)}
      
      if ( nrow(ae2) != 0 ){
        
        sub_ae2<- subset(ae2,i1==l) #l is end point
        #
        if (nrow(sub_ae2) ==0) {
          su_ae2<-NULL
          su_ae2<-as.data.frame(su_ae2)
          nrow(su_ae2)}
        #
        if (nrow(sub_ae2) >=1) {
          su_ae2<- subset(sub_ae2,targe==k) #k is start point
          
          if ( nrow(su_ae2) >=1 ) { sub_caus_ord<- subset(caus_ord, i==l) 
          su_caus_ord<- subset(sub_caus_ord, j==k) #compare with this (longest link)
          for(m in 1:nrow(su_ae2) ) {
            
            compare_caus_ord<- subset(caus_ord,j==su_ae2$targe[m])  # targe is start point
            compare_caus_or<- subset(compare_caus_ord,i==su_ae2$targe2[m])# i3 is end point, using this for next
            
            # compare with this
            if (compare_caus_or$max.rho_fromlag_realdata > su_caus_ord$max.rho_fromlag_realdata & compare_caus_or$lag > su_caus_ord$lag) {count_ae2<- count_ae2+1} } }
        }}
      
      count_all<- sum(count_ae+count_ae2)
      
      # merge
      #if (count_all >=1 ) {need_delete<-cbind(l,k) }   
      if ( count_all >=1 & count_all ==(nrow(su_ae)+nrow(su_ae2) ) ) {need_delete<-cbind(l,k) } 
      need_delet<-rbind(need_delet, need_delete )
      
    }
    
    need_deleted<-rbind(need_deleted, need_delet )
  }
  
  need_deleted<-as.data.frame(need_deleted)
  #unique
  if( nrow(need_deleted) >0 ){
    need_deleted<-unique(need_deleted) 
    colnames(need_deleted)<-c("end_point","start_point")
  }
  
  dim(need_deleted)
  write.csv(need_deleted,"need_deleted.csv")
  
  
  
  # test a number is in a vector
  # ae %in% ae2
  
  #=============== start remove indirect links
  # find it
  cau_ord<-read.csv("caus_ord_lag(mae).csv")
  which_row <- NULL
  if (nrow(need_deleted) > 0 ) {
    for (n in 1:dim(need_deleted)[1] ) {
      
      start_pointt <- need_deleted$start_point[n]
      end_pointt <- need_deleted$end_point[n]
      
      for (o in 1:dim(cau_ord)[1] ) { 
        if( cau_ord$j[o]==start_pointt & cau_ord$i[o]==end_pointt ) {which_row<- rbind(which_row,o)} }
      
    }
  }
  
  
  # remove it
  cau_or<-read.csv("caus_ord_lag(mae).csv")
  # 
  if (length(which_row) ==0) {cau_or2 <- cau_or}
  
  if (length(which_row) > 0) {cau_or2<-cau_or[-which_row,]} 
  
  
  
  for (ii in 1:length(cau_or2$i)) { cau_or2$total_links[ii]<-length(cau_or2$i[cau_or2$i==cau_or2$i[ii]])} # change totoal links
  dim(cau_or2) #211 10
  
  write.csv(cau_or2,"caus_ord_lag(mae)2.csv")
  
  # check midding species
  length( unique(cau_or2$i) ) #31, so missing one, the 18th species
  length( unique(cau_or2$j) ) #32
  
  
  
  
  # change 6
  ################################
  #### convert to ccm_sig; ccm_rho
  otago.links.data<-read.csv("caus_ord_lag(mae)2.csv")
  otago.graph.c<-graph.edgelist(as.matrix(otago.links.data[,c(4,5)]))
  ccm.sig<-get.adjacency(otago.graph.c,sparse=F)
  diag(ccm.sig)<-1
  
  ccm.rho<-matrix(0,nrow=number_species,ncol = number_species, byrow = T)
  otago.links.data<-as.data.frame( otago.links.data )
  for (k in unique(otago.links.data$i)) { 
    subs<-otago.links.data[otago.links.data$i==k,];
    ccm.rho[k,subs$j]<-subs$max.rho_fromlag_realdata
  }
  write.csv(ccm.sig,paste('ccm_sig_nin.csv',sep=''),row.names=F)
  write.csv(ccm.rho,paste('ccm_rho_nin.csv',sep=''),row.names=F)
  
  ccm.sig <- read.csv(paste('ccm_sig_nin.csv',sep=''),header=T,stringsAsFactors = F)
  ccm.rho <- read.csv(paste('ccm_rho_nin.csv',sep=''),header=T,stringsAsFactors = F)
  
  
  
  ###########  find ds ##############
  do <- read.csv('i.phyto.zoo.abu.biom.seas.csv',header=T,stringsAsFactors = F)
  
  Emax<-9
  out.sample <- F # T/F for out-of-sample forecast
  nout <- 0  # number of out-of-sample
  da.range <- 1:nrow(do) # Subsample for data analysis
  Ed<-E_alg_zoo[,2]
  
  dot <- do[da.range,1] # data time
  do <- do[da.range,-1] # time series of community data
  ndo <- nrow(do)
  nin <- ndo-nout # library sample size
  
  # Excluding rare species with no clear temporal dynamics 
  # end
  (nsp <- ncol(do))                            # number of selected species
  dim(do)
  # In-sample
  do.mean <- apply(do[1:nin,],2,mean,na.rm=T)  # mean abundance in in-sample
  do.sd <- apply(do[1:nin,],2,sd,na.rm=T)      # SD of abundance in in-sample
  d <- do[1:(nin-1),]                          # In-sample dataset at time t
  d_tp1 <- do[2:(nin),]                        # In-sample dataset at time t+1
  ds <- (d-repmat(do.mean,nrow(d),1))*repmat(do.sd,nrow(d),1)^-1 # Normalized in-sample dataset at time t
  ds_tp1 <- (d_tp1-repmat(do.mean,nrow(d_tp1),1))*repmat(do.sd,nrow(d_tp1),1)^-1 # Normalized in-sample dataset at time t+1
  dim(ds_tp1)
  # Out-sample
  if(out.sample&nout!=0){
    d.test <- do[nin:(ndo-1),]                 # Out-of-sample dataset at time t 
    dt_tp1 <- do[(nin+1):ndo,]                 # Out-of-sample dataset at time t+1
    ds.test <- (d.test-repmat(do.mean,nrow(d.test),1))*repmat(do.sd,nrow(d.test),1)^-1 # Normalized out-of-sample dataset at time t
    dst_tp1 <- (dt_tp1-repmat(do.mean,nrow(dt_tp1),1))*repmat(do.sd,nrow(dt_tp1),1)^-1 # Normalized out-of-sample dataset at time t+1
  }else{d.test <- dt_tp1 <- ds.test <- NULL}
  
  # Compiled data at time t 
  ds.all <- rbind(ds,ds.test)
  dim(ds)
  
  
  ######################################################################
  # Perform multiview embedding analysis for each node
  # Warning: It is time consuming for running multview embedding for each nodes
  do.multiview <- SaveFile <- T #T/F for saving files
  set.seed(49563)
  
  if(do.multiview){
    esele_lag <- esim.lag.demo(ds,ccm.rho,ccm.sig,Ed,kmax=1000,kn=100,max_lag=3,Emax=Emax)
    # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
    #original
    #if(SaveFile){write.csv(esele_lag,paste('eseleLag_',da.name,'_nin',nin,'_demo_NEW.csv',sep=''),row.names=F)}
    #me
    if(SaveFile){write.csv(esele_lag,paste('eseleLag_nin.csv',sep=''),row.names=F)}
  }
  
  #original
  #esele <- read.csv(paste('eseleLag_',da.name,'_nin',nin,'_demo.csv',sep=''),header=T)
  #me
  esele <- read.csv(paste('eseleLag_nin.csv',sep=''),header=T)
  
  ####################################################
  ## The computation of multiview distance
  dmatrix.mv <- mvdist.demo(ds,ds.all,esele)
  dmatrix.train.mvx <- dmatrix.mv[['dmatrix.train.mvx']]
  dmatrix.test.mvx <- dmatrix.mv[['dmatrix.test.mvx']]
  
  
  
  
  
  ######## Leave-one-out cross-validation for finding the optimal parameters for MDR S-map analysis
  ### Warning: The cross-validation is the most time-consuming step in MDR S-map requiring massive computations and .  
  ### Thus, we recommend dividing job into smaller parts (sub.da>1)  or used parallel computation (parall=T, ncore>=1)
  ### The example showing below divided the parameter space into five parts and ran independently (sub.da=5).
  do.MDR.CV <- F
  ### The parameter cv.unit determines the precision of selected parameters and strongly influences computation time.
  ### In our cases, we used cv.unit=0.025 to obtain more precise estimations
  ### This parameter may be adjusted to 0.05 or even 0.1, depending on how sensitive the results to parameter precision. 
  
  #original
  cv.unit <- 0.025
  #me
  cv.unit <- 0.3
  #original
  alpha.so <- seq(0.3, 1, cv.unit);            # Sequence of alpha
  
  sub.da <- 5                                # Divide the computation job into five parts 
  afsp <- eqsplit(1:length(alpha.so),sub.da) # Divide the parameter space based on alpha parameter
  # original
  #alf <- 5                                 # Run CV in the first parameter subset 
  for (alf in 1:1) { #me
    
    # Cross-validation of MDR analysis    
    #original
    #if(do.MDR.CV){
    #me
    if(do.MDR.CV==FALSE){
      #original
      #alpha.s <- alpha.so[afsp[alf,1]:afsp[alf,2]] # Subset parameter pace
      #me
      alpha.s <- alpha.so                           # Subset parameter pace
      #original
      #cv.ind <- cv.MDR.demo(ds, ds_tp1, dmatrix.list=dmatrix.train.mvx,parall=F, ncore=1, keep_intra=T,alpha.seq=alpha.s)
      #me
      cv.ind <- cv.MDR.demo(ds,ds_tp1,dmatrix.list=dmatrix.train.mvx, parall=T,ncore=1,keep_intra=T,lambda.seq = sort(logspace(-3,0,3),decreasing = T),
                            tht = c(3, 6, 9),
                            alpha.seq=alpha.so)
      
      # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
      #original
      #if(SaveFile){write.csv(cv.ind,paste('./Output/',da.name,'_nin',nin,'_cvunit',cv.unit,'_alph',alpha.s[1]*100,'_cvout_Nmvx_Rallx.csv',sep=''),row.names=F)}
      #me
      if(SaveFile){write.csv(cv.ind,paste('_nin_cvunit_alph.csv',sep=''),row.names=F)}
    }
    # Repeat the CV under different parameter subsets by changing alf=2,3..,sub.da                                  
    
  } # for loop (not orignal)
  
  
  
  
  ################################################################################
  # Compiled the CV results tested under different parts of parameter space
  #alpha.so[afsp[1:nrow(afsp),1]]*100
  
  CompileCV=T
  if(CompileCV){
    #original
    #cv.ind <- NULL; for(alf in 1:nrow(afsp)){
    #me
    cv.ind <- NULL; 
    #original
    #cv.ind <- rbind(cv.ind,read.csv(paste(da.name,'_nin',nin,'_cvunit',cv.unit,'_alph',alpha.so[afsp[alf,1]]*100,'_cvout_Nmvx_Rallx_demo.csv',sep=''),header=T))
    #me
    cv.ind <- read.csv(paste('_nin_cvunit_alph.csv',sep=''),header=T)
    
    # Select the optimal parameter set with the minimal MSE
    paracv.demo <- secv.demo(cv.ind)
    #original
    #write.csv(paracv.demo,paste(da.name,'_nin',nin,'_cvunit',cv.unit,'_OptimalCV_Nmvx_Rallx_NEW.csv',sep=''),row.names = F)
    #me
    write.csv(paracv.demo,paste('_OptimalCV_Nmvx_Rallx_NEW.csv',sep=''),row.names = F)
  }
  
  
  
  
  ############################################################
  # Fitting MDR S-map based on the parameters selected by CV
  #me
  do.MDR <- F
  #original
  do.MDR <- T
  #original
  cv.unit <- 0.025
  #me
  cv.unit <- 0.3                          
  ptype <- 'aenet'                           # enet:elastic-net or msaenet: adaptive elastic-net
  
  # Select the optimal parameter set with the minimal MSE
  #original
  #paracv.demo <- read.csv(paste(da.name,'_nin',nin,'_cvunit',cv.unit,'_OptimalCV_Nmvx_Rallx.csv',sep='')) #"model1024_0_0_0_nin98_cvunit0.025_OptimalCV_Nmvx_Rallx.csv"
  #me
  paracv.demo <- read.csv("_OptimalCV_Nmvx_Rallx_NEW.csv")
  
  if(do.MDR){
    # Fitting the MDR S-map
    smap.demo <- MDRsmap.demo(paracv=paracv.demo,ptype=ptype,keep_intra=T,out.sample=F,
                              ds,ds_tp1,ds.test,dst_tp1,
                              dmatrix.list=dmatrix.train.mvx,
                              dmatrix.test.list=dmatrix.test.mvx)
    
    # Save forecast skills
    nr.out <- smap.demo[['nr.out']];
    # To avoid overwrite the original files, we save them with different names, 'XXX_NEW'.
    if(SaveFile){
      #original
      #write.csv(nr.out,paste(da.name,'_nin',nin,'_cvunit',cv.unit,'_',ptype,'_nrout_Nmvx_Rallx_demo_NEW.csv',sep=''),row.names=F)
      #me
      write.csv(nr.out,paste('_nrout_Nmvx_Rallx_demo_NEW.csv',sep=''),row.names=F)
      # Save interaction Jacobian matrices at all time points
      #original
      #write.csv(smap.demo[['jcof']],paste(da.name,'_nin',nin,'_cvunit',cv.unit,'_',ptype,'_jcof_Nmvx_Rallx_demo_NEW.csv',sep=''),row.names=F)
      #me
      write.csv(smap.demo[['jcof']],paste('_jcof_Nmvx_Rallx_demo_NEW.csv',sep=''),row.names=F)
    }
  }
  
  
  
  #change 7
  #####################
  ### Jacobian
  Jocos <- read.csv(paste('_jcof_Nmvx_Rallx_demo_NEW.csv',sep='')) 
  
  eigen=correlationss=variance_inter=mean_int_variance_acrossss=mean_int_variance_withinss=mean_interaction_preda_prey_abs_mean=NULL
  mean_interaction_preda_prey_sqrt=min_interaction=median_interaction=diagnal1=diagnal2=diagnal_mean=traces=NULL
  evo_eco_evo_eigens=eco_evo_eco_eigens=D_evo_eco_evo_eigens=A_eco_evo_eco_eigens=D_matrix_eigens=A_matrix_eigens=NULL
  i=1
  for ( i in 1:nrow(ds)) {
    subs<-subset(Jocos,time==i)
    Jaco<-as.matrix( subs[,5:ncol(subs)] )
    
    if (any(is.na(Jaco))==TRUE) {deM <-  NA} 
    if (any(is.na(Jaco))==FALSE) {
      eigs <- eigen(Jaco)[["values"]]
      mx <- which.max(abs(eigs))
      deM <- abs(eigs)[mx]
      eigen<-rbind(eigen,deM)}
    #
    trac<-tr(Jaco)
    traces<-rbind(traces,trac)
    # correlation
    correlations1<-Jaco
    correlations1 <- as.matrix(correlations1)
    lowers<-correlations1[lower.tri(correlations1, diag = FALSE)]
    lowers2<-lowers[lowers!=0]
    uppers<-correlations1[upper.tri(correlations1, diag = FALSE)]
    uppers2<-uppers[uppers!=0]
    if ( length(uppers2)>length(lowers2) ){
      if ( any(is.na(  c(lowers2,uppers2) ))==TRUE ) {correlations2 <-  NA} 
      if ( any(is.na(  c(lowers2,uppers2) ))==FALSE ) {
        if (length(uppers2) <3) {correlations2 <-  NA}
        if (length(uppers2) >=3) {correlations2 <- cor.test(uppers2, c(lowers2,rep(0,times=length(uppers2)-length(lowers2))), method = "pearson")$estimate } }}
    if ( length(uppers2)<length(lowers2) ){
      if ( any(is.na(  c(lowers2,uppers2) ))==TRUE ) {correlations2 <-  NA} 
      if ( any(is.na(  c(lowers2,uppers2) ))==FALSE ) {
        if (length(lowers2) <3) {correlations2 <-  NA}
        if (length(lowers2) >=3) {correlations2 <- cor.test( c(uppers2,rep(0,times=length(lowers2)-length(uppers2))), lowers2, method = "pearson")$estimate } }}
    if ( length(uppers2)==length(lowers2) ){
      if ( any(is.na(  c(lowers2,uppers2) ))==TRUE ) {correlations2 <-  NA} 
      if ( any(is.na(  c(lowers2,uppers2) ))==FALSE ) {correlations2 <- cor.test(uppers2, lowers2, method = "pearson")$estimate } }
    correlationss<-rbind(correlationss,correlations2)
    
    # variance in predator.prey 
    mean_int_variance<-Jaco
    mean_int_variance <- as.matrix(mean_int_variance)
    diag( mean_int_variance ) <- 0  
    variance_inte<-var(mean_int_variance[1:(nrow(Jaco)*nrow(Jaco))])
    variance_inter<-rbind(variance_inter,variance_inte)
    
    # predator.prey 
    mean_int<-Jaco
    mean_int <- as.matrix(mean_int)
    diag( mean_int ) <- 0  
    mean_interactio_preda_prey_abs_mean<-sum( mean_int[mean_int!=0] ) / length(mean_int[mean_int!=0])
    if (any(is.na(mean_interactio_preda_prey_abs_mean))==TRUE) {mean_interactio_preda_prey_abs_mean <-NA} 
    mean_interaction_preda_prey_abs_mean<-rbind(mean_interaction_preda_prey_abs_mean,mean_interactio_preda_prey_abs_mean)
    
    
    #diagna
    mean_diag<-Jaco
    mean_diag <- as.matrix(mean_diag)
    diagnal_mea<- sum( diag(mean_diag) )/ length( diag(mean_diag) ) # just mean
    diagnal_mean<-rbind(diagnal_mean,diagnal_mea)
    
    #
    A_matrix<- matrix( Jaco[1:23,1:23], nrow=23, byrow=F)
    B_matrix<- matrix( Jaco[1:23,24:46], nrow=23, byrow=F)
    C_matrix<- matrix( Jaco[24:46,1:23], nrow=23, byrow=F)
    D_matrix<- matrix( Jaco[24:46,24:46], nrow=23, byrow=F)
    
    evo_eco_evo<- C_matrix %*% inv(A_matrix) %*% (-B_matrix)
    eco_evo_eco<- B_matrix %*% inv(D_matrix) %*% (-C_matrix)  
    D_evo_eco_evo<-D_matrix+evo_eco_evo
    A_eco_evo_eco<-A_matrix+eco_evo_eco
    
    #evo_eco_evo
    if (any(is.na(evo_eco_evo))==TRUE) {evo_eco_evo_eigen <-  NA} 
    if (any(is.na(evo_eco_evo))==FALSE) {
      eigs <- eigen(evo_eco_evo)[["values"]]
      mx <- which.max(abs(eigs))
      evo_eco_evo_eigen <- abs(eigs)[mx]
      evo_eco_evo_eigens<-rbind(evo_eco_evo_eigens,evo_eco_evo_eigen) }
    #D_evo_eco_evo
    if (any(is.na(D_evo_eco_evo))==TRUE) {D_evo_eco_evo_eigen <-  NA} 
    if (any(is.na(D_evo_eco_evo))==FALSE) {
      eigs <- eigen(D_evo_eco_evo)[["values"]]
      mx <- which.max(abs(eigs))
      D_evo_eco_evo_eigen <- abs(eigs)[mx]
      D_evo_eco_evo_eigens<-rbind(D_evo_eco_evo_eigens,D_evo_eco_evo_eigen) }
    #A_eco_evo_eco
    if (any(is.na(A_eco_evo_eco))==TRUE) {A_eco_evo_eco_eigen <-  NA} 
    if (any(is.na(A_eco_evo_eco))==FALSE) {
      eigs <- eigen(A_eco_evo_eco)[["values"]]
      mx <- which.max(abs(eigs))
      A_eco_evo_eco_eigen <-  abs(eigs)[mx]
      A_eco_evo_eco_eigens<-rbind(A_eco_evo_eco_eigens,A_eco_evo_eco_eigen) }
    #eco_evo_eco
    if (any(is.na(eco_evo_eco))==TRUE) {eco_evo_eco_eigen <-  NA} 
    if (any(is.na(eco_evo_eco))==FALSE) {
      eigs <- eigen(eco_evo_eco)[["values"]]
      mx <- which.max(abs(eigs))
      eco_evo_eco_eigen <- abs(eigs)[mx]
      eco_evo_eco_eigens<-rbind(eco_evo_eco_eigens,eco_evo_eco_eigen) }   
    # A matrix
    if (any(is.na(A_matrix))==TRUE) {A_matrix_eigen <-  NA} 
    if (any(is.na(A_matrix))==FALSE) {
      eigs <- eigen(A_matrix)[["values"]]
      mx <- which.max(abs(eigs))
      A_matrix_eigen <- abs(eigs)[mx]
      A_matrix_eigens<-rbind(A_matrix_eigens,A_matrix_eigen) }
    # D matrix
    if (any(is.na(D_matrix))==TRUE) {D_matrix_eigen <-  NA} 
    if (any(is.na(D_matrix))==FALSE) {
      eigs <- eigen(D_matrix)[["values"]]
      mx <- which.max(abs(eigs))
      D_matrix_eigen <- abs(eigs)[mx]
      D_matrix_eigens<-rbind(D_matrix_eigens,D_matrix_eigen) }
    
  } #loop  1:nrow(ds)
  
  eigen<-as.data.frame(cbind(eigen,evo_eco_evo_eigens,D_evo_eco_evo_eigens,eco_evo_eco_eigens,A_matrix_eigens,D_matrix_eigens,A_eco_evo_eco_eigens))
  colnames(eigen)<-c("eigen","evo_eco_evo","D_evo_eco","eco_evo_eco","A_matrix_eigens","D_matrix_eigens","A_eco_evo")
  head(eigen)
  
  plot(eigen[,1],type = "l", ylim=c(0.8,20))
  abline(h=1, col="purple")
  eigen[,1][which(eigen[,1]<1)]   
  length(eigen[,1][which(eigen[,1]<1)])
  length(eigen[,1][which(eigen[,1]>1)])
  #colnames(eigen[,1])<-c("eigenvalue")
  write.csv(eigen,"eigen.csv")
  
  
  #
  plot(traces,type = "l", ylim=c(0.8,5))
  colnames(traces)<-c("traces")
  write.csv(traces,"traces.csv")
  ### 
  plot(correlationss,type = "l", ylim=c(-1,1))
  abline(h=0)
  write.csv(correlationss,"correlation_preda_prey_no_abs.csv")
  
  # vairance in predator-prey
  plot(variance_inter,type = "l", ylim=c(-0.05,0.05))
  abline(h=0)
  write.csv(variance_inter,"variance_preda_prey_no_abs.csv")
  
  #preda_prey_abs_mean
  plot(mean_interaction_preda_prey_abs_mean,type = "l", ylim=c(0,0.2))
  abline(h=0)
  write.csv(mean_interaction_preda_prey_abs_mean,"mean_interaction_preda_prey_abs_mean.csv")
  #
  ###  diag 3
  plot(diagnal_mean,type = "l", ylim=c(-1,1))
  write.csv(diagnal_mean,"diagnal_geome_mean.csv")
  
  #
  mean_interaction_all<-cbind(correlationss,variance_inter,
                              mean_interaction_preda_prey_abs_mean,diagnal_mean)
  colnames(mean_interaction_all)<-c("correlation_predator_prey",
                                    "variance_predator_prey",
                                    "preda_prey_abs_mean","diagnal_geomean")
  write.csv(mean_interaction_all,"mean_interaction_all.csv")
  
  
  
  
  
  
  #### computing diversity index
  new_alg_zo<-read.csv("i.phyto.zoo.abu.biom.seas.csv")
  new_alg_zo<-as.data.frame(new_alg_zo)# having year
  dim(new_alg_zo)
  # scale
  alga_zo<-new_alg_zo[,2:dim(new_alg_zo)[[2]]] # no year
  dim(alga_zo) # 161  27
  head(alga_zo) # 161 19
  colnames(alga_zo)

  #
  library(vegan)
  
  # phyto+zoo
  #### richness
  species_richness<-NULL
  for (i in 1:dim(alga_zo)[[1]]) {species_richness[i] <- sum(alga_zo[i,]>0) }
  species_richness
  # shannon, simpson
  shan_all<-diversity(alga_zo, index = "shannon") 
  sim_all<-diversity(alga_zo, index = "simpson") 
  
  shan_sim_all<- cbind(shan_all,sim_all,species_richness)
  colnames(shan_sim_all)<-c("shan_al","sim_al","species_richness")
  write.csv(shan_sim_all,"shan_sim_only.csv")
  
  #---- combine eigen,interaction-diversity-----
  #----#----#----#----#----#----#----#----#----
  New_year<-read.csv("i.physical_Shackelton Point.csv")
  day_temp<-cbind(as.data.frame(New_year$sampledate), New_year$Avg.Temp)
  colnames(day_temp)<-c("date","temp_station_all")

  shan_sim_only_final<-cbind(day_temp, shan_sim_all, rbind(eigen,NA), rbind(traces,NA), rbind(mean_interaction_all,NA))
  write.csv(shan_sim_only_final,"shan_sim_all2.csv")
  head(shan_sim_only_final)
  
  
  
  
  
  
  
  ######################
  ####lyapunov_exponent
  Jocos <- read.csv(paste('_jcof_Nmvx_Rallx_demo_NEW.csv',sep='')) 
  head(Jocos)
  lyapunov_exponents=NULL
  i=1
  for ( i in 1:5) {  #for ( i in 1:6) {
#  for ( i in 1:11) {  
    sub<-subset(Jocos,time==i)
    subs<-sub[5:ncol(sub)]
    
    #whole matrix
    whole_matri<-as.matrix( subs[,1:ncol(subs)] )
    if (i==1){whole_matrix<- diag(ncol(whole_matri))}
    #whole_matri[whole_matri==0]<- 0.008
    whole_matrix<-whole_matrix %*% whole_matri
    dim(whole_matrix)
    dim(whole_matri)
    
    #Ecology matrix 
    Ecology_matri<-as.matrix( subs[1:(nrow(subs)/2),1:(ncol(subs)/2)] )
    if (i==1){Ecology_matrix<- diag(ncol(Ecology_matri))}
    #Ecology_matri[Ecology_matri==0]<- 0.008#
    Ecology_matrix<-Ecology_matrix %*% Ecology_matri
    dim(Ecology_matrix)
    dim(Ecology_matri)
    
    #Evolution matrix 
    Evolution_matri<-as.matrix( as.matrix( subs[(nrow(subs)/2+1):(nrow(subs)), (ncol(subs)/2+1):(ncol(subs))] ) )
    if (i==1){Evolution_matrix<- diag(ncol(Evolution_matri))}
    #Evolution_matri[Evolution_matri==0]<- 0.008#
    Evolution_matrix<-Evolution_matrix %*% Evolution_matri
    dim(Evolution_matrix)
    dim(Evolution_matri)
    
    #A-B(D-1)C
    B_matrix<-as.matrix( subs[1:(nrow(subs)/2),(ncol(subs)/2+1):ncol(subs)] )
    C_matrix<-as.matrix( subs[(ncol(subs)/2+1):ncol(subs),1:(nrow(subs)/2)] )
    A.plus.eco.evo.eco<-Ecology_matri-B_matrix %*% inv(Evolution_matri)  %*% C_matrix
    if (i==1){A.plus.eco.evo.eco_matrix<- diag(ncol(A.plus.eco.evo.eco))}
    A.plus.eco.evo.eco_matrix<-A.plus.eco.evo.eco_matrix %*% A.plus.eco.evo.eco
    
    #D+C(A-1)B
    B_matrix<-as.matrix( subs[1:(nrow(subs)/2),(ncol(subs)/2+1):ncol(subs)] )
    C_matrix<-as.matrix( subs[(ncol(subs)/2+1):ncol(subs),1:(nrow(subs)/2)] )
    D.plus.evo.eco.evo<-Evolution_matri-C_matrix %*% inv(Ecology_matri)  %*% B_matrix
    if (i==1){D.plus.evo.eco.evo_matrix<- diag(ncol(D.plus.evo.eco.evo))}
    D.plus.evo.eco.evo_matrix<-D.plus.evo.eco.evo_matrix %*% D.plus.evo.eco.evo
  } 
  
  
  # whole max(eigen)
  eigs_whole <- eigen(whole_matrix)[["values"]]
  mx_whole <- which.max(abs(eigs_whole))
  deM_whole <- abs(eigs_whole)[mx_whole]
  
  # ecology max(eigen)
  eigs_Ecology <- eigen(Ecology_matrix)[["values"]]
  mx_Ecology <- which.max(abs(eigs_Ecology))
  deM_Ecology <- abs(eigs_Ecology)[mx_Ecology]
  
  # Evolution max(eigen)
  eigs_Evolution <- eigen(Evolution_matrix)[["values"]]
  mx_Evolution <- which.max(abs(eigs_Evolution))
  deM_Evolution <- abs(eigs_Evolution)[mx_Evolution]
  
  #A-B(D-1)C max(eigen)
  eigs_A.plus.eco.evo.eco <- eigen(A.plus.eco.evo.eco_matrix)[["values"]]
  mx_A.plus.eco.evo.eco <- which.max(abs(eigs_A.plus.eco.evo.eco))
  deM_A.plus.eco.evo.eco <- abs(eigs_A.plus.eco.evo.eco)[mx_A.plus.eco.evo.eco]
  
  #D+C(A-1)B max(eigen)
  eigs_D.plus.evo.eco.evo <- eigen(D.plus.evo.eco.evo_matrix)[["values"]]
  mx_D.plus.evo.eco.evo <- which.max(abs(eigs_D.plus.evo.eco.evo))
  deM_D.plus.evo.eco.evo <- abs(eigs_D.plus.evo.eco.evo)[mx_D.plus.evo.eco.evo]
  
  
  #lyapunov exponents
  t<-max(Jocos$time)
  LE_whole<- log( abs(det(deM_whole * as.matrix(whole_matrix))) ) /t
  LE_whole
  LE_Ecology<-log( abs(det(deM_Ecology * as.matrix(Ecology_matrix))) ) /t
  LE_Ecology
  LE_Evolution<-log( abs(det(deM_Evolution * as.matrix(Evolution_matrix))) ) /t
  LE_Evolution
  
  #LE[A-B(D-1)C]
  LE_A.plus.eco.evo.eco<-LE_whole+log(deM_A.plus.eco.evo.eco/(deM_whole*abs(det(Evolution_matrix)))) /t
  LE_A.plus.eco.evo.eco
  #LE[D+C(A-1)B]  
  LE_D.plus.evo.eco.evo<-LE_whole+log(deM_D.plus.evo.eco.evo/(deM_whole*abs(det(Ecology_matrix)))) /t
  LE_D.plus.evo.eco.evo
  
  #
  shan_sim_only_final$LE_whole<-rep(LE_whole,times=nrow(shan_sim_only_final))
  shan_sim_only_final$LE_Ecology<-rep(LE_Ecology,times=nrow(shan_sim_only_final))
  shan_sim_only_final$LE_Evolution<-rep(LE_Evolution,times=nrow(shan_sim_only_final))
  
  #
  shan_sim_only_final$LE_A.plus.eco.evo.eco<-rep(LE_A.plus.eco.evo.eco,times=nrow(shan_sim_only_final))
  shan_sim_only_final$LE_D.plus.evo.eco.evo<-rep(LE_D.plus.evo.eco.evo,times=nrow(shan_sim_only_final))
  
  write.csv(shan_sim_only_final,"shan_sim_all2.csv")
  head(shan_sim_only_final)
  
  
  
  
  
  
  