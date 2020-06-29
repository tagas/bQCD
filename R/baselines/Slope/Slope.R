#!/usr/bin/Rscript



##### Infer the best model

## Slope: Working title for algorithm

##Params:
# a + bx + cx^2 + dx^3 + e*exp(x) + fx^-1 + gx^-2 + h*log2(x)
fofx = function(x, p){
    x.neg = x
    x.neg[x == 0] = resolution
    y.head = p$a + p$b * x + p$c * (x^2) + p$d * (x^3) + p$e * (exp(x)) + p$f * (x.neg^(-1)) + p$g * (x.neg^(-2))
    return(y.head)
}
getFunctionIndex = function(p){
    if(abs(p$b) > 0){
        return(1)
    }else if(abs(p$c) > 0){
        return(2)
    }else if(abs(p$d) > 0){
        return(3)
    }else if(abs(p$e) > 0){
        return(4)
    }else if(abs(p$f) > 0){
        return(5)
    }else if(abs(p$g) > 0){
        return(6)
    }else{
        return(0)
    }
}
###### Best fit
fitLin = function(y,x){
    m = lm(y~x)
    params = list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
    sse = tail(anova(m)[,2],1)
    coeff = naTo0(m$coefficients)
    params$a = coeff[1]
    params$b = coeff[2]
    res = list(sse=sse, model=parameterScore(m), par=params)
    return(res)
}
fitQuad = function(y,x){
    m = lm(y~I(x^2))
    params = list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
    sse = tail(anova(m)[,2],1)
    coeff = naTo0(m$coefficients)
    params$a = coeff[1]
    params$c = coeff[2]
    res = list(sse=sse, model=parameterScore(m), par=params)
    return(res)
}
fitCub = function(y,x){
    m = lm(y~I(x^3))
    params = list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
    sse = tail(anova(m)[,2],1)
    coeff = naTo0(m$coefficients)
    params$a = coeff[1]
    params$d = coeff[2]
    res = list(sse=sse, model=parameterScore(m), par=params)
    return(res)
}
fitExp = function(y,x){
    xe = exp(x)
    m = lm(y~xe)
    params = list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
    sse = tail(anova(m)[,2],1)
    coeff = naTo0(m$coefficients)
    params$a = coeff[1]
    params$e = coeff[2]
    res = list(sse=sse, model=parameterScore(m), par=params)
    return(res)
}
fitNegLin = function(y,x){
    x[x == 0] = resolution
    xe = x^(-1)
    m = lm(y~xe)
    params = list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
    sse = tail(anova(m)[,2],1)
    coeff = naTo0(m$coefficients)
    params$a = coeff[1]
    params$f = coeff[2]
    res = list(sse=sse, model=parameterScore(m), par=params)
    return(res)
}
fitNegQuad = function(y,x){
    x[x == 0] = resolution
    xe = x^(-1)
    m = lm(y~xe)
    params = list(a=0, b=0, c=0, d=0, e=0, f=0, g=0, h=0)
    sse = tail(anova(m)[,2],1)
    coeff = naTo0(m$coefficients)
    params$a = coeff[1]
    params$g = coeff[2]
    res = list(sse=sse, model=parameterScore(m), par=params)
    return(res)
}

findBestFit = function(y,x){
    r1 = fitLin(y,x)
    r2 = fitQuad(y,x)
    r3 = fitCub(y,x)
    r4 = fitExp(y,x)
    r5 = fitNegLin(y,x)
    r6 = fitNegQuad(y,x)
    s = rep(0,6)
    s[1] = gaussian_score_emp_sse(r1$sse, length(x)) + r1$model
    s[2] = gaussian_score_emp_sse(r2$sse, length(x)) + r2$model
    s[3] = gaussian_score_emp_sse(r3$sse, length(x)) + r3$model
    s[4] = gaussian_score_emp_sse(r4$sse, length(x)) + r4$model
    s[5] = gaussian_score_emp_sse(r5$sse, length(x)) + r5$model
    s[6] = gaussian_score_emp_sse(r6$sse, length(x)) + r6$model
    i = which(s == min(s))
    if(i[1] == 1){
        return(r1)
    }else if(i[1] == 2){
        return(r2)
    }else if(i[1] == 3){
        return(r3)
    }else if(i[1] == 4){
        return(r4)
    }else if(i[1] == 5){
        return(r5)
    }else{
        return(r6)
    }
}

fitWrapperStirr = function(y,x){
    maxX = max(x)
    xf = round(x*100000)
    tx = table(xf)
    if(length(tx) <= 2){ ## binary data
        return(list(s=0, b=-2, p=-2))
    }
    mx = mean(tx)
    sd = sd(tx)
    costsND = 2147483647
    potential_bins = 0
    if(mx >= 10){
        potential_bins = 1
        xi = xf %in% as.numeric(names(tx)[tx >= 5])
        xfg = xf[xi]
        xg = x[xi]
        yg = y[xi]
        score1 = 0
        score2 = 0
        score3 = 0
        score4 = 0
        score5 = 0
        score6 = 0
        sse = 0.0
        if(sum(!xi) > 1){
            xl = x[!xi]
            yl = y[!xi]
            resSingletons = findBestFit(yl,xl)
            sse = gaussian_score_emp_sse(resSingletons$sse, length(xl)) + resSingletons$model
        }
        for(e in unique(sort(xfg))){
            ones = xfg == e
            yt = sort(yg[ones])
            xt = normX(1:length(yt), 10) - 5.0
            f1 = fitLin(yt, xt)
            f2 = fitQuad(yt, xt)
            f3 = fitCub(yt, xt)
            f4 = fitExp(yt, xt)
            f5 = fitNegLin(yt, xt)
            f6 = fitNegQuad(yt, xt)
            score1 = score1 + gaussian_score_emp_sse(f1$sse, length(xt)) + f1$model
            score2 = score2 + gaussian_score_emp_sse(f2$sse, length(xt)) + f2$model
            score3 = score3 + gaussian_score_emp_sse(f3$sse, length(xt)) + f3$model
            score4 = score4 + gaussian_score_emp_sse(f4$sse, length(xt)) + f4$model
            score5 = score5 + gaussian_score_emp_sse(f5$sse, length(xt)) + f5$model
            score6 = score6 + gaussian_score_emp_sse(f6$sse, length(xt)) + f6$model
        }
        score = sse + min(score1,score2,score3,score4,score5,score6)
        modelCosts = 1 + log2nChoosek(length(x)-1, length(unique(sort(xfg)))-1) + logg(length(unique(sort(xfg)))) ## cut points
        costs = score + modelCosts
        return(list(s=costs, b=c(1)))
        if(costs < costsND){
            costsND = costs
        }
    }
    res = findBestFit(y,x)
    score = gaussian_score_emp_sse(res$sse, length(x)) + res$model
    modelCosts = 1 ## additional model cost of 1 cause no bins have to be identified
    costsD = score + modelCosts
    if(costsND < costsD){
        return(list(s=costsND, b=1, p=potential_bins))
    }else{
        return(list(s=costsD, b=0, p=potential_bins))
    }
}

fitComparison = function(fit, l, lx, score, bins_old){
    newScore = gaussian_score_emp_sse(fit$sse, l) + fit$model
    delta = newScore - score + log2nChoosek(lx-1, bins_old+1) + logg(bins_old+1) - logg(bins_old) - log2nChoosek(lx-1, bins_old)
    return(delta)
}

fitWrapper = function(y,x){
    minNum = 5
    maxX = max(x)
    xf = round(x*100000)
    tx = table(xf)
    if(length(tx) <= 2){  ## binary data
        return(list(s=0, b=-2, p=-2))
    }
    mx = mean(tx)
    sd = sd(tx)
    lx = length(x)
    res = findBestFit(y,x)
    fun = getFunctionIndex(res$par)
    score = gaussian_score_emp_sse(res$sse, length(x))
    modelCosts = 1 + res$model ## additional model cost of 1 cause no bins have to be identified
    costs = score + modelCosts
    xi = xf %in% as.numeric(names(tx)[tx >= minNum])
    xfg = xf[xi]
    xg = x[xi]
    yg = y[xi]
    score1 = score
    score2 = score
    score3 = score
    score4 = score
    score5 = score
    score6 = score
    bins = rep(1,6)
    for(e in unique(sort(xfg))){
        ones = xfg == e
        yt = sort(yg[ones])
        xt = normX(1:length(yt), 10) - 5.0
        curr_l = length(xt)
        old_sse = sum((yg[ones] - fofx(xg[ones], res$par))^2)
        old_score = gaussian_score_emp_sse(old_sse, curr_l)
        f1 = fitLin(yt, xt)
        f2 = fitQuad(yt, xt)
        f3 = fitCub(yt, xt)
        f4 = fitExp(yt, xt)
        f5 = fitNegLin(yt, xt)
        f6 = fitNegQuad(yt, xt)
        s1 = fitComparison(f1, curr_l, lx, old_score, bins[1])
        if(s1 < 0){
            bins[1] = bins[1] + 1
            score1 = score1 + s1
        }
        s2 = fitComparison(f2, curr_l, lx, old_score, bins[2])
        if(s2 < 0){
            bins[2] = bins[2] + 1
            score2 = score2 + s2
        }
        s3 = fitComparison(f3, curr_l, lx, old_score, bins[3])
        if(s3 < 0){
            bins[3] = bins[3] + 1
            score3 = score3 + s3
        }
        s4 = fitComparison(f4, curr_l, lx, old_score, bins[4])
        if(s4 < 0){
            bins[4] = bins[4] + 1
            score4 = score4 + s4
        }
        s5 = fitComparison(f5, curr_l, lx, old_score, bins[5])
        if(s5 < 0){
            bins[5] = bins[5] + 1
            score5 = score5 + s5
        }
        s6 = fitComparison(f6, curr_l, lx, old_score, bins[6])
        if(s6 < 0){
            bins[6] = bins[6] + 1
            score6 = score6 + s6
        }
    }
    # correct for all divided
    for(i in 1:length(bins)){
        if(bins[i] > length(tx)){
            bins[i] = bins[i] - 1
        }
    }
    # add model costs
    score1 = score1 + 1 + log2nChoosek(length(x)-1, bins[1]-1)
    score2 = score2 + 1 + log2nChoosek(length(x)-1, bins[2]-1)
    score3 = score3 + 1 + log2nChoosek(length(x)-1, bins[3]-1)
    score4 = score4 + 1 + log2nChoosek(length(x)-1, bins[4]-1)
    score5 = score5 + 1 + log2nChoosek(length(x)-1, bins[5]-1)
    score6 = score6 + 1 + log2nChoosek(length(x)-1, bins[6]-1)
    scores = c(score1,score2,score3,score4,score5,score6)
    i = which(scores == min(scores))[1]
    score = scores[i]
    costs = score
    potential_bins = sum(tx >= minNum)
    if(potential_bins != length(tx)){
        potential_bins = potential_bins + 1
    }
    if(bins[i] > 1){
        fun = i
    }
    return(list(s=score, b=bins[i], p=potential_bins, f=fun))
}

fitFG = function(y,x){
    xf = round(x*100000)
    tx = table(xf)
    mx = mean(tx)
    ## if non-deterministic
    if(mx >= 10){
        res = findBestFit(y,x)
        model = res$model
        y.rem = y - fofx(x, res$par)
        xt = normX(1:length(y), 10) - 5.0
        res2 = findBestFit(sort(y.rem),xt)
        model = model + res$model
        data = gaussian_score_emp_sse(res2$sse, length(x))
        return(data + model)
    }else{
        return(fitWrapper(y,x))
    }
}

defaultScore = function(x){
    score = -logg(resolution) * length(x)
    return(score)
}

### Slope algorithm
Slope = function(t, prune=1.0, alpha=0.001){
    ## remove percentage of outliers
    if(prune < 1.0){
        require("ldbod")
        count = floor(dim(t)[1] * prune)
        scores = ldbod(t, k=10)
        t <- t[order(scores$lof,decreasing=F)[1:count],]
    }
    x = normX(t[,1],1)
    y = normX(t[,2],1)
    print("Calculate X->Y...")
    setResolution(mindiff(y))
    print(resolution)
    dy = defaultScore(y)
    resXtoY = fitWrapper(y,x)
    sseXtoY = resXtoY$s
    
    print("Calculate Y->X...")
    setResolution(mindiff(x))
    print(resolution)
    dx = defaultScore(x)
    resYtoX = fitWrapper(x,y)
    sseYtoX = resYtoX$s
    
    dXY = sseXtoY + dx
    dYX = sseYtoX + dy
    dXtoY = dXY / (dx + dy)
    dYtoX = dYX / (dx + dy)
    
    # Get delta
    eps = dXtoY - dYtoX
    pv = 2^(-(abs(dXY - dYX)/2))
    ## if data is binary
    if(resXtoY$p == -2 | resYtoX$p == -2){
        eps = 0.0
        pv = 1.0
    }
    
    # Determine causal direction
    causd = "--"
    if(abs(eps) > 0.0 & pv < alpha){
        if(eps < 0){
            causd = "->"
        }else{
            causd = "<-"
        }
    }

    r = list(epsilon = eps, cd = causd, p.value=pv)
    return(r)
}

SlopeInfo = function(t, prune=1.0, alpha=0.001){
    ## remove percentage of outliers
    if(prune < 1.0){
        require("ldbod")
        count = floor(dim(t)[1] * prune)
        scores = ldbod(t, k=10)
        t <- t[order(scores$lof,decreasing=F)[1:count],]
    }
    x = normX(t[,1],1)
    y = normX(t[,2],1)
    #print("Calculate X->Y...")
    setResolution(mindiff(y))
    #print(resolution)
    dy = defaultScore(y)
    resXtoY = fitWrapper(y,x)
    sseXtoY = resXtoY$s
    
    #print("Calculate Y->X...")
    setResolution(mindiff(x))
    #print(resolution)
    dx = defaultScore(x)
    resYtoX = fitWrapper(x,y)
    sseYtoX = resYtoX$s
    
    dXY = sseXtoY + dx
    dYX = sseYtoX + dy
    dXtoY = dXY / (dx + dy)
    dYtoX = dYX / (dx + dy)
    
    # Get delta
    eps = dXtoY - dYtoX
    pv = 2^(-(abs(dXY - dYX)/2))
    ## if data is binary
    if(resXtoY$p == -2 | resYtoX$p == -2){
        eps = 0.0
        pv = 1.0
    }
    
    # Determine causal direction
    causd = "--"
    if(abs(eps) > 0.0 & pv < alpha){
        if(eps < 0){
            causd = "->"
        }else{
            causd = "<-"
        }
    }
    
    r = list(epsilon = eps, cd = causd, p.value=pv, sc=c(dXY, dYX))
    return(r)
}
