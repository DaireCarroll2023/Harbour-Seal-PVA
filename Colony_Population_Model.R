#Daire Carroll University of Gothenburg 2022
#Grey seal PVA script, functions for predictive modelling

#################################
#Packages
library(ggplot2)
library(data.table)
#################################

#t = number of cycles
#A = maximum age
#Kc = vector of carrying capacities 
#P = initial population as an age structured matrix
#B = age specific birthrates as a vector
#S = age specific survival as a vector
#Rc = random variation terms as a vector
#Hc = size of hunt as a vector, either proportion or set quota
#H_delay = Numnber of cycles before hunting is implemented, can be 0... numeric
#H_sk = hunt skew, a vector expressing how the hunt should be applied to the population
#proportion = TRUE/FALSE if hunt should be applied as a proportion of the population (TRUE) or as a set quota (FALSE) 

#P_D = Probability of a disease outbreak in each cycle... numeric
#O_D = Outcome of disease, maximum proportion of an age class killed during an outbreak... numeric
#W_D = Proportion of O_D experienced by age class... vector

structure = function(B,S,N){
  B = B*S/2
  L = matrix(B,ncol = A)
  for(i in 1:(A-1)){
    x = numeric(A)
    x[i] = S[i]
    L = rbind(L, x)
  }
  ev = eigen(L)
  ev = as.numeric(ev$vectors[,1])
  P = c()
  for(i in 1:A){
    P[i] = ev[i]/sum(ev)*N
  } 
  return(P)
} #find dominant left eigen vector for given demographic scenario, with birthrates B (vector), survival S (vector) and popualiton size

eig_val = function(B,S){
  A = length(B)
  B = B*S/2
  L = matrix(B,ncol = A)
  for(i in 1:(A-1)){
    x = numeric(A)
    x[i] = S[i]
    L = rbind(L, x)
  }
  ev = eigen(L)
  return(as.numeric(ev$value[1]))
} #find dominant left eigen value for given demographic scenario

PopMod = function(t, A, Kc, P, B, S, R, Hc, H_delay, H_sk, proportion, P_D, O_D, W_D){
  
  ###housekeeping and troubleshooting begin###
  
  if(t==0){		
    
    if(stru == TRUE){
      
      return(P)
      
    }else{
      
      return(sum(P))
      
    }
    
  } 
  
  null_args = list(H_sk,E_sk)
  null_args_names = c("H_sk","E_sk")
  for(i in 1:length(null_args)){
    if(is.null(null_args[[i]])){
      assign(paste0(null_args_names[i]), 0)
    }
  }
  
  if(nrow(P)!=A || (length(S))!=A || length(B)!=A){
    stop('Lengths of P, S, & B must = A')
  }
  
  ###housekeeping and troubleshooting end###
  
  ###function for Leslie matrix assembly begin###
  
  const_leslie = function(B,S){
    
    B = B*S/2
    
    L = matrix(B,ncol = A) 
    
    for(i in 1:(A-1)){
      x = numeric(A)
      x[i] = S[i]
      L = rbind(L, x)
    }	
    
    return(L)
    
  }
  
  L = const_leslie(B,S)
  ev = eigen(L)
  ev1 = as.numeric(ev$value[1]) #dominant eigan value (lambda), also 'long term population growth'
  
  ###function for Leslie matrix assembly end###
  
  ###carrying capacity begin### 
  
  Keffect = function(K,ev1,TP){ 
    
    if(is.null(K)){
      return(1)
    }else{
      return(((K+(ev1-1)*TP)/K)^-1)
    }  
  } #calculate the effect of carrying capactity at t based on population size at t-1, currently for best L1 matrix dominent eigen value
  
  ###carrying capacity end###
  
  ###hunting function start###
  
  even_spread = function(PI){
    remaining = sqrt(sum(PI[which(PI < 0)])^2)
    PI[which(PI < 0)] = 0               
    PI = PI - remaining*(PI/sum(PI))
    return(PI)
  }  
  
  if(proportion == TRUE){
    
    apply_H = function(PI,H,H_sk){
      
      basic_H = function(PI,H){
        
        PI = PI-PI*H
        return(PI)
        
      } #in cases where there is no skew hunting is applied evenly across age classes
      
      if(sum(H_sk) == 0){ 
        PI = basic_H(PI,H)
        return(PI) #with no hunting skew, basic hunting is applied across classes
      }else{
        
        PI = PI = PI-sum(PI)*H*H_sk
        PI = even_spread(PI)
        
        return(PI) #in cases where hunting is applied as a vector, when an age class is depleted, the quota is met with individuals from other classes
        
      }
      
    }
    
  }else if(proportion == FALSE){
    
    apply_H = function(PI,H,H_sk){
      
      if(sum(PI)<H){
        PI[1:length(PI)] = 0
        return(PI)
      }
      
      basic_H = function(PI,H){
        
        PI = PI-H*(PI/sum(PI))
        return(PI)
        
      } #in cases where there is no skew hunting is applied evenly across age classes
      
      PI = even_spread(PI)
      
      if(sum(H_sk) == 0){ 
        
        PI = basic_H(PI,H)
        
        return(PI) #with no hunting skew, basic hunting is applied across classes
        
      }else{
        
        PI = PI - H*H_sk
        PI = even_spread(PI)
        
        return(PI) #in cases where hunting is applied as a vector, when an age class is depleted, the quota is met with individuals from other classes
        
      }
      
    }
  }
  
  check = function(P){
    P[P < 0] = 0
    return(P)
  }
  
  ###hunting function end###
  
  ###stochastic function begin###
  
  stocst = function(R,L){
    stoc = (1 - (rbeta(1,5,5)*2-1)*R)*L
    return(stoc)    
  } #random variation applied to L
  
  ###stochastic function end###
  
  ###Migration function begin####
  
  migration = function(P_out,M){
    
    P_mod = matrix(0, nrow = nrow(P_out),ncol = ncol(P_out))
    
    for(i in 1:ncol(P_out)){
      for(j in 1:ncol(P_mod)){
        P_mod[,j] = P_mod[,j] + P_out[,i]*M[i,j]
      }
    }  
    
    return(P_mod)
    
  } #migration function
  
  ###Migration function end####

  ####Holding matrices begin####  

  P_fin = matrix(0, nrow = (t+1),ncol = ncol(P)) #this will be the final output - a matrix with total population size for each of the colonies for each t
  P_fin[1,] = colSums(P)
  names(P_fin) = names(P) 
  P_out = matrix(0, nrow = nrow(P),ncol = ncol(P)) #a holding matrix for the structured popualiotns at t
  
  ####Holding matrices end####
  
  ###first year begin###

  for(i in 1:(ncol(P))){
    
    P0 = P[,i]
    K = Kc[i]
    H = Hc[i]
    
    TP = sum(P0)
    Kef = Keffect(K,ev1,TP)
    LK = L*Kef #apply effect of carrying capacity to Leslie matrix
    LK = stocst(R,LK)
    PI = LK%*%P0	 #population at t = 1 
    
    if(H > 0 && H_delay < 1){ 
      PI = apply_H(PI,H,H_sk) #effects of hunting
    } #apply hunting if delay less than 1 cycle is specified
    
    PI = check(PI)
    
    P_out[,i] = PI
    
  }

  if(sum(M)>0){
    P_out = migration(P_out,M)
  } #apply migration
  
  if(P_D > 0 && runif(1, 0, 1)<P_D){
    impact = W_D*O_D
    P_out = P_out - P_out*impact
    check(P_out)
  } #apply disease
  
  P_fin[2,] = colSums(P_out)
  
  if(t==1){
    return(P_fin) #end the function if only a single year is being run
  }
  
  H_delay = H_delay - 1
  
  ###first year end###

  ####cycle through remaining years###
    
  for(j in 2:t){
      
    for(i in 1:(ncol(P))){
        
        P0 = P_out[,i]
        K = Kc[i]
        H = Hc[i]
        
        TP = sum(P0)
        Kef = Keffect(K,ev1,TP)
        LK = L*Kef #apply effect of carrying capacity to Leslie matrix
        LK = stocst(R,LK)
        PI = LK%*%P0	 #population at t = 1 
        
        if(H > 0 && H_delay < 1){
          PI = apply_H(PI,H,H_sk) #effects of hunting
        } #apply hunting if delay less than 1 cycle is specified
        
        PI = check(PI)
        
        P_out[,i] = PI
        
      }
      
    if(sum(M)>0){
      P_out = migration(P_out,M)
    } #apply migration
    
    if(P_D > 0 && runif(1, 0, 1)<P_D){
      impact = W_D*O_D
      P_out = P_out - P_out*impact
      check(P_out)
    } #apply disease
    
    P_fin[(j + 1),] = colSums(P_out)
    
    H_delay = H_delay - 1
      
    }
    
  return(P_fin)
  
} #Population model with colony structure, lets think about the order in which things occor

#################################
#standard parameters

t = 50 
A = 38

Koster = c(0,281,223,178,141,120,102,86,74,62,54,45,38,33,28,24,20,17,14,13,11,9,8,7,6,4,4,4,3,2,3,1,2,1,1,1,1,1)
Vaderoarna = c(0,162,129,102,81,69,59,50,43,36,31,26,22,19,16,14,12,10,9,7,6,5,5,4,4,3,3,3,2,2,2,1,2,1,1,1,1,1)
Lysekil = c(0,137,109,86,69,58,50,42,36,30,26,22,19,16,14,12,10,9,7,6,5,5,4,4,3,2,2,2,2,1,2,1,1,1,1,1,1,1)
Marstrand = c(0,260,207,164,131,111,95,80,69,57,50,42,35,31,26,22,19,16,13,12,10,8,7,7,6,4,4,4,3,2,3,1,2,1,1,1,1,1)
Onsala = c(0,296,235,187,149,126,108,91,78,65,57,48,40,35,29,25,21,18,15,13,11,9,8,7,6,4,4,4,3,2,3,1,2,1,1,1,1,1)
Varberg = c(0,161,128,102,81,69,59,50,43,36,31,26,22,19,16,14,12,10,9,7,6,5,5,4,4,3,3,3,2,2,2,1,2,1,1,1,1,1)
Hallvadero = c(0,156,124,99,78,67,57,48,41,35,30,25,21,19,15,13,11,10,8,7,6,5,5,4,4,3,3,3,2,2,2,1,2,1,1,1,1,1)
Laso = c(0,106,84,67,53,45,39,33,28,24,20,17,15,13,11,9,8,7,6,5,4,4,3,3,3,2,2,2,2,1,2,1,1,1,1,1,1,1)
Hesselo = c(0,194,154,123,98,83,71,60,51,43,37,31,26,23,19,17,14,12,10,9,8,6,6,5,4,3,3,3,2,2,2,1,2,1,1,1,1,1)
Anholt = c(0,107,85,68,54,46,39,33,29,24,21,18,15,13,11,9,8,7,6,5,4,4,3,3,3,2,2,2,2,1,2,1,1,1,1,1,1,1)
SWKattegat = c(0,126,100,79,63,54,46,39,33,28,24,20,17,15,13,11,9,8,7,6,5,4,4,3,3,2,2,2,2,1,2,1,1,1,1,1,1,1)  #initial population

P =  cbind(Koster,Vaderoarna,Lysekil,Marstrand,Onsala,Varberg,Hallvadero,Laso,Hesselo,Anholt,SWKattegat)

B = rep(0,A)
B[4] = 0.17*2
B[5] = 0.33*2
B[6:27] = 0.47*2
B[27:(A-1)] = 0.35*2

S = rep(0,A)
S[1] = 0.75
S[2:4] = 0.89
S[5:(A-1)] = 0.95 

R = 0.05 

P_D = 0.01
O_D = 0.65 
W_D = c(100,11.55,9.95,9.15,9,10.25,9,8.5,7.25,7,5.6,2.6,1.9,0.2,3,1.7,0.1,0,0.2,1.1,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11)/100 #weighting of disease per age class

M_Koster=c(98,1,1,0,0,0,0,0,0,0,0) 
M_Vaderoarna=c(1,97,1,1,0,0,0,0,0,0,0)
M_Lysekil=c(1,1,96,1,1,0,0,0,0,0,0)
M_Marstrand=c(0,1,1,95,1,1,1,0,0,0,0)
M_Onsala=c(0,1,1,1,95,1,1,0,0,0,0)
M_Varberg=c(0,0,0,1,1,95,1,1,1,0,0)
M_Hallvadero=c(0,0,0,1,1,1,96,1,1,1,1)
M_Laso=c(0,0,0,0,0,1,1,96,1,1,1)
M_Hesselo=c(0,0,0,0,0,0,1,1,96,1,1)
M_Anholt=c(0,1,0,0,0,0,1,1,1,95,1)
M_SWKattegat=c(0,0,0,0,0,0,1,1,1,1,96) 

M =rbind(M_Koster,M_Vaderoarna,M_Lysekil,M_Marstrand,M_Onsala,M_Varberg,M_Hallvadero,M_Laso,M_Hesselo,M_Anholt,M_SWKattegat)/100

colnames(M) = colnames(P)
rownames(M) = colnames(P)

H_sk = NULL

H_delay = 0

####################################################