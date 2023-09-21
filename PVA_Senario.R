#################################
#Daire Carroll, University of Gothenburg, 2021, daire.carroll@bioenv.gu.se
#################################
#Packages
if(!require("ggplot2")){
  install.packages("ggplot2")
}
if(!require("reshape2")){
  install.packages("reshape2")
}

library(ggplot2)
library(reshape2)

#################################
#parameters for harbour seals, Silva et al. 2021

Ka=c(2380,1370,1155,2204,2509,1363,1319,894,1645,905,1062) 

t = 100 

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
B[4] = 0.17
B[5] = 0.33
B[6:27] = 0.47
B[27:(A-1)] = 0.35

S = rep(0,A)
S[1] = 0.75
S[2:4] = 0.89
S[5:(A-1)] = 0.95 

f = 0.1  

R = 0.05 

H = 0.025 
H_delay = 5 

Target_ages = NULL 
H_skew = NULL  
H_d = NULL

P_D = 0.07
O_D = 0.65 
W_D = c(100,11.55,9.95,9.15,9,10.25,9,8.5,7.25,7,5.6,2.6,1.9,0.2,3,1.7,0.1,0,0.2,1.1,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11)/100 #weighting of disease per age class

E = 0

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

#################################
#run model 'reps' number of times saving the mean curve and quasi extinction risk

reps = 100 #number of iterations, takes roughly 
QE = 500 #quasi extinction value

t_tot_vals = matrix(0, ncol = ncol(P), nrow = (t+1))
Prob_QE = rep(0, ncol(P))

ptm = proc.time()

for(i in 1:reps){
  
  tot_vals = c()
  
  for(j in 1:ncol(P)){
    tot_vals[j] = (sum(P[,j]))
  }
  
  tot_vals = matrix(tot_vals, ncol = ncol(P))
  
  colnames(tot_vals) = colnames(P)
  
  for(j in 1:t){
    
    Pt = PopMod(t = j,
                A = A,
                Ka = Ka,
                P = P,
                B = B,
                S = S,
                H = H,
                f = f,
                H_delay =  H_delay,
                R = R,
                Target_ages = NULL,
                H_skew = NULL,
                H_d = NULL,
                P_D = P_D,
                O_D = O_D,
                W_D = W_D,
                M = M,
                E = E
                )
    
    Tot_Pop = c()
    
    for(j in 1:ncol(Pt)){
      Tot_Pop[j] = sum(Pt[,j])
    }
    
    tot_vals = rbind(tot_vals,Tot_Pop)
    
  }
  
  for(j in 1:ncol(tot_vals)){
    if(length(which(tot_vals[,j] < QE))>1){
      Prob_QE[j] = Prob_QE[j]+1
    } #if the population dips below QE at any point it is counted as extinction... can be replaced with poulations state after t: if(length(tot_vals[,j][length(tot_vals[,j])] < QE))>1){Prob_QE[j] = Prob_QE[j]+1}
  }  
  
  for(j in 1:ncol(P)){
    
    t_tot_vals[,j] = t_tot_vals[,j] + tot_vals[,j]
    
  }
  
}

Prob_QE = Prob_QE/reps
Prob_QE = rbind(Prob_QE)
colnames(Prob_QE) = colnames(P)

m_tot_vals = t_tot_vals/reps
colnames(m_tot_vals) = colnames(P)

proc.time() - ptm

#################################
#Plot average output

cleanup_grid = theme(panel.border = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(color = "black"),
)

cleanup_text = theme(axis.text.x = element_text(color = "black", size = 12, angle = 0, hjust = .5, vjust = .5, face = "plain"),
                     axis.text.y = element_text(color = "black", size = 12, angle = 0, hjust = 0.8, vjust = 0.4, face = "plain"),  
                     axis.title.x = element_text(color = "black", size = 20, angle = 0, hjust = .5, vjust = 0, face = "plain"),
                     axis.title.y = element_text(color = "black", size = 20, angle = 90, hjust = .5, vjust = .5, face = "plain"))		

m_tot_vals2 = melt(m_tot_vals, varnames = c("Year","Population"))

p1 = ggplot(m_tot_vals2, aes(x = Year, y = value, colour = Population)) +
  geom_line() +
  cleanup_grid +
  cleanup_text +
  labs(x="Year", y=expression("Population Size"))+
  geom_hline(yintercept = QE, col = "red", linetype = "dashed")

p1
print(Prob_QE)
