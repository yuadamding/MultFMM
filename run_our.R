{
  setwd("E:/OneDrive - Binghamton University/Simulations")
  library("MASS")
  library("mvtnorm")
  source("helper.R")
  source("fma.init_3.R")
  source("mixFMM_07_26.R")
  source("functions.R")
  data = read.csv(
    'C:\\Users\\Adam\\Dropbox\\meeting_minutes_Yu_Ding\\HCHS\\data_na_omited_07_11_2022.csv'
  )
  selected = ifelse(as.numeric(data[, colnames(data) == 'AGE']) >= 18 &
                      as.numeric(data[, colnames(data) == 'AGE']) <= 24, 1, 0)
  data_selected = data[selected == 1, ]
  col_selected = c(
    'PID' ,
    'AGE' ,
    'RACE' ,
    'GENDER' ,
    'MARITAL_STATUS' ,
    'INCOME' ,
    'EMPLOYED' ,
    'SLEA12A' ,
    'SLEA12B' ,
    'SLEA12C' ,
    'SLEA12D' ,
    'SLEA12E' ,
    'SLEA12F' ,
    'SLEA12G' ,
    'SLEA12H' ,
    'SLEA12I' ,
    'SLEA12J' ,
    'SLEA4' ,
    'SLEA5' ,
    'SLEA6' ,
    'SLEA7' ,
    'SLEA11' ,
    'SLPA54' ,
    'SLPA60' ,
    'SLPA90' ,
    'SLPA91' ,
    'SLPA92' ,
    'SLPA93' ,
    'SLPA103' ,
    'SLPA111' ,
    'SLPA112' ,
    'SLPA113' ,
    'SLPA114' ,
    'SLPA115' ,
    'SLPA116' ,
    'SCEA1' ,
    'SCEA2' ,
    'SCEA3' ,
    'SCEA4' ,
    'SCEA5' ,
    'SCEA6' ,
    'SCEA7' ,
    'SCEA8' ,
    'SCEA9' ,
    'SCEA10' ,
    'WBEA1' ,
    'WBEA2' ,
    'WBEA3' ,
    'WBEA4' ,
    'WBEA5' ,
    'WBEA6' ,
    'WBEA7' ,
    'WBEA8' ,
    'WBEA9' ,
    'WBEA10' ,
    'WBEA11' ,
    'WBEA12' ,
    'WBEA13' ,
    'WBEA14' ,
    'WBEA15' ,
    'WBEA16' ,
    'WBEA17' ,
    'WBEA18' ,
    'WBEA19' ,
    'WBEA20' ,
    'GPAQ_WORK_VIG' ,
    'GPAQ_WORK_MOD' ,
    'GPAQ_TRSPORT' ,
    'GPAQ_REC_VIG' ,
    'GPAQ_REC_MOD',
    'AHEI1',
    'AHEI2',
    'AHEI3',
    'AHEI4',
    'AHEI5',
    'AHEI6',
    'AHEI7',
    'AHEI8',
    'AHEI9',
    'AHEI10',
    'AHEI11' ,
    'LABA66' ,
    'LABA68' ,
    'LABA70' ,
    'LABA71' ,
    'LABA72' ,
    'LABA74' ,
    'LABA75' ,
    'LABA76' ,
    'LABA91' ,
    'LABA101' ,
    'LABA103' ,
    'LABA67' ,
    'LABA69' ,
    'INSULIN_FAST' ,
    'HOMA_IR' ,
    'HBA1C_SI' ,
    'BMI' ,
    'WAIST_HIP' ,
    'ALCOHOL_USE' ,
    'CIGARETTE_USE' ,
    'SBPA5' ,
    'SBPA6'
  )
  data_selected$GENDER = ifelse(data_selected$GENDER == 'F' , 1, 0)
  X1 = as.matrix(data_selected[, 3:8])
  Xc1 = list(as.matrix(data_selected[, 9:18]))
  Xc2 = list(as.matrix(data_selected[, 19:23]))
  X2 = (as.matrix(data_selected[, 24:36]))
  Xc3 = list(as.matrix(data_selected[, 37:46]))
  Xc4 = list(as.matrix(data_selected[, 47:56]))
  Xc5 = list(as.matrix(data_selected[, 57:66]))
  X3 = (as.matrix(data_selected[, 67:71]))
  X4 = (as.matrix(data_selected[, 72:82]))
  X5 = (as.matrix(data_selected[, 83:104]))
  p_list = array(c(6, 13, 5, 11, 22), 5)
  X = array(0, c(863, 22, 4))
  standardization = function(x) {
    shape = dim(x)
    for (i in 1:shape[2]) {
      sigma = sd(x[, i])
      mu = mean(x[, i])
      for (j in 1:shape[1]) {
        x[j, i] = (x[j, i] - mu) / sigma
      }
    }
    return(x)
  }
  #X[, 1:p_list[1], 1] = standardization(X1)
  X[, 1:p_list[2], 1] = standardization(X2)
  X[, 1:p_list[3], 2] = standardization(X3)
  X[, 1:p_list[4], 3] = standardization(X4)
  X[, 1:p_list[5], 4] = standardization(X5)
  Xc = c(Xc1, Xc2, Xc3, Xc4, Xc5)
  #Xc = c(Xc1, Xc2)
}
load("1_exp6_Xc.RData")
Xc1 = Xc
load("2_exp6_Xc.RData")
Xc2 = Xc
load("3_exp6_Xc.RData")
Xc3 = Xc
load("4_exp6_Xc.RData")
Xc4 = Xc
load("5_exp6_Xc.RData")
Xc5 = Xc
load("X_6.RData")
Xc = c(Xc1, Xc2, Xc3, Xc4, Xc5)
res = mixfmm_07_26(
  X = X,
  Xc = Xc,
  k = 5,
  lambda_cnlm = 20,
  lambda_cnlm2 =20
)


j_list = array(0, 40)
for (i in 1:dim(Xc1[[1]])[2]) {
  j_list[i] = max(Xc1[[1]][, i])
}
X = array(0, c(400, 40, 2))
X[,  , 1] = Xc1[[1]]
X[,  , 2] = Xc1[[1]]
Xc = c(Xc1, Xc1)
res = mixfmm_07_26(
  X = NULL,
  Xc = Xc,
  k = 3,
  lambda_cnlm = 20,
  lambda_cnlm2 = 20
)

#---------------------------------------------------
library("doSNOW")
library("foreach")
lambda_cnlm_seq = seq(1, 101, 7)
lambda_cnlm2_seq = seq(1, 101, 7)

mc = 5
cl = makeCluster(11)
registerDoSNOW(cl)

bic = array(dim = c(length(lambda_cnlm_seq), length(lambda_cnlm2_seq)))
time = array(dim = c(length(lambda_cnlm_seq), length(lambda_cnlm2_seq)))
muc = array(dim = c(length(lambda_cnlm_seq), length(lambda_cnlm2_seq), 3, mc * 2))
beta = array(dim = c(length(lambda_cnlm_seq), length(lambda_cnlm2_seq), 3, 2, 5))

findIndx = function(i, len1, len2) {
  num_2 = ifelse(i %% len2 == 0, len2, i %% len2)
  num_1 = i %/% len2 + ifelse(i %% len2 == 0, 0, 1)
  return(c(num_1, num_2))
}

print(Sys.time())
output = foreach(i = 1:(length(lambda_cnlm_seq) * length(lambda_cnlm2_seq)), .export = "mixfmm_07_26") %dopar%
  {
    num_1 = findIndx(i, length(lambda_cnlm_seq), length(lambda_cnlm2_seq))[1]
    num_2 = findIndx(i, length(lambda_cnlm_seq), length(lambda_cnlm2_seq))[2]
    lambda_cnlm = lambda_cnlm_seq[num_1]
    lambda_cnlm2 = lambda_cnlm2_seq[num_2]
    fit = mixfmm_07_26(
      X = X,
      Xc = Xc,
      k = 3,
      lambda_cnlm = lambda_cnlm,
      lambda_cnlm2 = lambda_cnlm2
    )
    list(num_1, num_2, fit)
  }
print(Sys.time())
save(output, file = "res_exp_6.RData")
cat(" Finish! \n")
stopCluster(cl)

for (i in 1:(length(lambda_cnlm_seq) * length(lambda_cnlm2_seq))) {
  num_1 = output[[i]][[1]]
  num_2 = output[[i]][[2]]
  bic[num_1, num_2] = output[[i]][[3]]$bic
  time[num_1, num_2] = output[[i]][[3]]$time
  muc[num_1, num_2, ,] = output[[i]][[3]]$muc
  beta[num_1, num_2, , ,] = output[[i]][[3]]$beta
}
stopCluster(cl)

save(bic, file = "exp7_bic.RData")
save(time, file = "exp7_time.RData")
save(muc, file = "exp7_muc.RData")
save(beta, file = "exp7_beta.RData")
save(output, file = "exp7_output.RData")
# stopCluster(cl)
#
require(plot3D)
persp3D(
  x = lambda_cnlm_seq,
  y = lambda_cnlm2_seq,
  z = bic,
  theta = 225,
  phi = 25,
  expand = 0.75,
  xlab = "lambda_cnlm_seq",
  ylab = "lambda_cnlm2_seq",
  ticktype = "detailed",
  zlab = "BIC",
  axes = TRUE
)


plot(bic[, 2])
