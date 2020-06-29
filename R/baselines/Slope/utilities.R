#!/usr/bin/Rscript

sysouts = T



### na to zero
naTo0 = function(c){
    c[is.na(c)] = 0
    return(c)
}

###### Normalization
normX = function(x, n){
    if(min(x) == max(x)){
        return( rep(n, length(x)) )
    }else{
        return( ((x-min(x)) / (max(x) - min(x))) * n )
    }
}

logg = function(x){
    if(x == 0){
        return(0)
    }else{
        return(log2(x))
    }
}
log2fac = function(n){
    sum = 0
    for(i in 2:n){
        sum = sum + logg(i)
    }
    return(sum)
}
log2nChoosek = function(n, k){
    if(k > n | k == 0){
        return(0)
    }else{
        return(log2fac(n) - log2fac(k) - log2fac(n-k))
    }
}
logN = function(z){
    z = ceiling(z)
    if(z < 1){
        return(0)
    }else{
        logstar = logg(z)
        sum = logstar
        while(logstar > 0){
            logstar = logg(logstar)
            sum = sum + logstar
        }
        return(sum + logg(2.865064))
    }
}

fitTestWMW = function(y,x,cycles=100,alpha=0.05){
    set.seed(1234)
    sigCount = 0
    for(i in 1:cycles){
        random_permutation = sample(length(x))
        cp = ceiling(length(x) / 2)
        x.train = x[random_permutation[1:cp]]
        y.train = y[random_permutation[1:cp]]
        x.test  = x[random_permutation[(cp+1):length(x)]]
        y.test  = y[random_permutation[(cp+1):length(x)]]
        res = findBestFit(y.train,x.train)
        err.train = (y.train - fofx(x.train, res$par))^2
        err.test  = (y.test  - fofx(x.test,  res$par))^2
        wmw.res = wilcox.test(err.train, err.test, alternative="two.sided")
        p.value = wmw.res$p.value
        if(!is.na(p.value) & p.value < alpha){
            sigCount = sigCount + 1
        }
    }
    return(sigCount/cycles)
}

getDuplicatePositions = function(x){
    last = x[1]
    pos = 0
    s = ""
    for(i in 2:length(x)){
        if(x[i] != last){
            last = x[i]
            if(i - (pos + 1) > 1){
                k = paste(pos, (i-1), sep=":")
                if(nchar(s) > 0){
                    s = paste(s, k, sep=";")
                }else{
                    s = k
                }
            }
            pos = i - 1
        }
    }
    # last elem involved
    if((pos + 1) != length(x)){
        k = paste(pos, length(x), sep=":")
        if(nchar(s) > 0){
            s = paste(s, k, sep=";")
        }else{
            s = k
        }
    }
    return(s)
}

###### Entropy and conditional entropy
entropy = function(x){
    x.unique = unique(sort(x))
    N = length(x)
    H = 0
    for(i in 1:length(x.unique)){
        frac = sum(x == x.unique[i]) / N
        H = H - (frac * logg(frac))
    }
    return(H)
}
conditionalEntropy = function(x,y){
    y.unique = unique(sort(y))
    N = length(y)
    H = 0
    for(i in 1:length(y.unique)){
        xiy = x[y == y.unique[i]]
        H = H + ( (length(xiy) / N) * entropy(xiy) )
    }
    return(H)
}

###### Discretization
equiWidthBinning = function(dataset, bins){
    newd = cut(dataset, bins)
    return(as.numeric(newd))
}

###### Helper
getMeanOccurance = function(){
    df = data.frame(I=uv, X1=rep(0,length(uv)), SD1=rep(0,length(uv)), X2=rep(0,length(uv)), SD2=rep(0,length(uv)), R=ref.uv$V6)
    for(i in 1:length(uv)){
        t = readI(uv[i])
        df$X1[i] = mean(table(t[,1]))
        df$SD1[i] = sd(table(t[,1]))
        df$X2[i] = mean(table(t[,2]))
        df$SD2[i] = sd(table(t[,2]))
    }
    return(df)
}
getMinOccurance = function(){
    df = data.frame(I=uv, X1=rep(0,length(uv)), SD1=rep(0,length(uv)), X2=rep(0,length(uv)), SD2=rep(0,length(uv)), R=ref.uv$V6)
    for(i in 1:length(uv)){
        t = readI(uv[i])
        df$X1[i] = min(table(t[,1]))
        df$SD1[i] = sd(table(t[,1]))
        df$X2[i] = min(table(t[,2]))
        df$SD2[i] = sd(table(t[,2]))
    }
    return(df)
}
getDims = function(){
    df = data.frame(I=uv, L=rep(0,length(uv)))
    for(i in 1:length(uv)){
        t = readI(uv[i])
        df$L[i] = dim(t)[1]
    }
    return(df)
}

###### Calculate decision rate (-- = 1/2)
decision_rate_w = function(res){
    return(decision_rate(res$Correct, res$Eps, res$Cds))
}
decision_rate = function(corr, eps, cds){
    df = data.frame(A=corr, B=abs(eps), C=cds)
    df = df[with(df, order(-B)),]
    sum = 0
    dr = rep(0, dim(df)[1])
    for(i in 1:dim(df)[1]){
        if(df$C[i] == "--"){
            sum = sum + 0.5
        }else{
            sum = sum + df$A[i]
        }
        dr[i] = sum / i
    }
    return(c(1,dr))
}
decision_rate_meta = function(res, uv=uv, oneHalf=T){
    corr = res$Correct
    eps = res$Eps
    cds = res$Cds
    df = data.frame(Uv=uv, A=corr, B=abs(eps), C=cds)
    df = df[with(df, order(-B)),]
    sum = 0
    total = 0
    step = rep(0, dim(df)[1])
    dr = rep(0, dim(df)[1])
    for(i in 1:dim(df)[1]){
        val = meta$V6[df$Uv[i]]
        total = total + val
        if(df$C[i] == "--"){
            sum = sum + 0.5 * val
            if(oneHalf){
                dr[i] = sum / total
            }else{
                dr[i] = dr[i-1]
            }
        }else{
            sum = sum + df$A[i] * val
            dr[i] = sum / total
        }
        step[i] = total
    }
    step = step / total
    return(list(D=c(1,dr), S=c(0,step)))
}
decision_rate_meta_pos = function(res, uv=uv){
    corr = res$Correct
    eps = res$Eps
    cds = res$Cds
    df = data.frame(Uv=uv, A=corr, B=abs(eps), C=cds)
    df = df[with(df, order(B)),]
    sum = 0
    total = 0
    step = rep(0, dim(df)[1])
    dr = rep(0, dim(df)[1])
    for(i in 1:dim(df)[1]){
        val = meta$V6[df$Uv[i]]
        total = total + val
        if(df$C[i] == "--"){
            sum = sum + 0.5 * val
        }else{
            sum = sum + df$A[i] * val
        }
        dr[i] = sum / total
        step[i] = total
    }
    step = step / total
    return(list(D=c(1,dr), S=c(0,step)))
}
decision_rate_w_pos = function(res){
    corr = res$Correct
    eps = res$Eps
    cds = res$Cds
    df = data.frame(A=corr, B=abs(eps), C=cds)
    df = df[with(df, order(B)),]
    sum = 0
    dr = rep(0, dim(df)[1])
    for(i in 1:dim(df)[1]){
        if(df$C[i] == "--"){
            sum = sum + 0.5
        }else{
            sum = sum + df$A[i]
        }
        dr[i] = sum / i
    }
    return(c(1,dr))
}

###### SSE
SSE = function(x, yhead){
    e = sum((x-yhead)^2)
    return(e)
}

###### PlotI
plotI = function(i){plot(readI(i))}

###### MDL score
resolution = 0.01
setResolution = function(val){
    resolution <<- val
}
gaussian_score_emp = function(x){
    sse = SSE(x, mean(x))
    var = sse / length(x)
    sigma = sqrt(var)
    return(gaussian_score(sigma, x))
}
gaussian_score = function(sigma, x){
    sse = SSE(x, mean(x))
    n = length(x)
    sigmasq = sigma^2
    if(sse == 0.0 | sigmasq == 0.0){
        return(0.0)
    }else{
        err = (sse / (2 * sigmasq * log(2))) + ((n/2) * logg(2 * pi * sigmasq)) - n * logg(resolution)
        return(err)
    }
}
gaussian_score_emp_sse = function(sse, n){
    var = sse / n
    sigma = sqrt(var)
    return(gaussian_score_sse(sigma, sse, n))
}

gaussian_score_sse = function(sigma, sse, n){
    sigmasq = sigma^2
    if(sse == 0.0 | sigmasq == 0.0){
        return(0.0)
    }else{
        err = (sse / (2 * sigmasq * log(2))) + ((n/2) * logg(2 * pi * sigmasq)) - n * logg(resolution)
        return(max(err,0))
    }
}
parameterScore = function(model){
    sum = 0
    coeff = naTo0(model$coefficients)
    for(c in coeff){
        if(c != 0){
            ca = abs(c)
            cdummy = ca
            prec = 1
            while(cdummy < 1000){
                cdummy = cdummy * 10
                prec = prec + 1
            }
            sum = sum + logN(cdummy) + logN(prec) + 1
        }
    }
    return(sum)
}

### igci
S = function(x.uns){
    sum = 0.0
    x = sort(x.uns)
    m = length(x)
    for(i in 1:(m-1)){
        sum = sum + logg(abs(x[i+1] - x[i]))
    }
    sum = sum * (1/(m-1))
    sum = sum + logg(m)
    return(sum)
}
normalize = function(tt){
    t = tt
    # normalize to mean 0 and sd 1
    for(j in 1:(dim(t)[2])){
        t[,j] = as.numeric(t[,j])
        t[,j] = (t[,j] - min(t[,j])) / (max(t[,j]) - min(t[,j]))
    }
    return(t)
}


normalize_pobs = function(tt){
  n <- nrow(tt)
  u1 <- rank(as.numeric(tt[,1]), ties.method = "random")/(n + 1)
  u2 <- rank(as.numeric(tt[,2]),ties.method = "random")/(n + 1)
  return(cbind(u1,u2))
}



# case 1
# % uniform reference measure
# x = (x - min(x)) / (max(x) - min(x));
# y = (y - min(y)) / (max(y) - min(y));
# case 2
# % Gaussian reference measure
# x = (x - mean(x)) ./ std(x);
# y = (y - mean(y)) ./ std(y);
# otherwise

normalize_G = function(tt){
  t = tt
  # normalize to mean 0 and sd 1
  for(j in 1:(dim(t)[2])){
    t[,j] = as.numeric(t[,j])
    #t[,j] = (t[,j] - min(t[,j])) / (max(t[,j]) - min(t[,j]))
    t[,j] = (t[,j] - mean(t[,j]) / sd(t[,j]))
  }
  return(t)
}

## mindiff
mindiff = function(x){
    xs = sort(x)
    diff = 0.01
    for(i in 1:(length(x)-1)){
        new_diff = xs[i+1] - xs[i]
        #print(new_diff)
        #print(xs[i+1])
        if(new_diff != 0 & new_diff < diff){
            diff = new_diff
        }
    }
    return(diff)
}
data_precision = function(x){
    precision = 1
    set = x != round(x)
    x = x[set]
    while(length(x) > 0){
        precision = precision / 10
        x = 10 * x
        set = x != round(x)
        x = x[set]
    }
    return(precision)
}

##### IGCI
IGCI = function(t){
    t = normalize(t)#normalize(t)
    diffs = S(t[,2]) - S(t[,1])
    causd = "--"
    if(diffs < 0){
        causd = "->"
    }else if(diffs > 0){
        causd = "<-"
    }
    return(list(cd=causd, epsilon=diffs))
}

##### IGCI
IGCI_G = function(t){
  t = normalize_G(t)
  diffs = S(t[,2]) - S(t[,1])
  causd = "--"
  if(diffs < 0){
    causd = "->"
  }else if(diffs > 0){
    causd = "<-"
  }
  return(list(cd=causd, epsilon=diffs))
}

