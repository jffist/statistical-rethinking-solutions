library(rethinking)

#the only way to fix "couldn't coerse S4 object to double" error while plotting rethink objects
#functions were extracted from
#https://github.com/rmcelreath/rethinking/search?utf8=%E2%9C%93&q=setmethod%28%22plot&type=Code


setMethod("plot" , "compareIC" , function(x,
                                          y,
                                          xlim,
                                          SE = TRUE,
                                          dSE = TRUE,
                                          weights = FALSE,
                                          ...) {
  dev_in <- x@output[[1]] - x@output[[2]] * 2
  dev_out <- x@output[[1]]
  if (!is.null(x@output[['SE']]))
    devSE <- x@output[['SE']]
  dev_out_lower <- dev_out - devSE
  dev_out_upper <- dev_out + devSE
  if (weights == TRUE) {
    dev_in <- ICweights(dev_in)
    dev_out <- ICweights(dev_out)
    dev_out_lower <- ICweights(dev_out_lower)
    dev_out_upper <- ICweights(dev_out_upper)
  }
  n <- length(dev_in)
  if (missing(xlim)) {
    xlim <- c(min(dev_in), max(dev_out))
    if (SE == TRUE & !is.null(x@output[['SE']])) {
      xlim <- c(min(dev_in), max(dev_out_upper))
    }
  }
  main <- colnames(x@output)[1]
  set_nice_margins()
  dotchart(
    dev_in[n:1] ,
    labels = rownames(x@output)[n:1] ,
    xlab = "deviance" ,
    pch = 16 ,
    xlim = xlim ,
    ...
  )
  points(dev_out[n:1] , 1:n)
  mtext(main)
  # standard errors
  if (!is.null(x@output[['SE']]) & SE == TRUE) {
    for (i in 1:n) {
      lines(c(dev_out_lower[i], dev_out_upper[i]) , rep(n + 1 - i, 2) , lwd =
              0.75)
    }
  }
  if (!all(is.na(x@dSE)) & dSE == TRUE) {
    # plot differences and stderr of differences
    dcol <- col.alpha("black", 0.5)
    abline(v = dev_out[1] , lwd = 0.5 , col = dcol)
    diff_dev_lower <- dev_out - x@output$dSE
    diff_dev_upper <- dev_out + x@output$dSE
    if (weights == TRUE) {
      diff_dev_lower <- ICweights(diff_dev_lower)
      diff_dev_upper <- ICweights(diff_dev_upper)
    }
    for (i in 2:n) {
      points(
        dev_out[i] ,
        n + 2 - i - 0.5 ,
        cex = 0.5 ,
        pch = 2 ,
        col = dcol
      )
      lines(
        c(diff_dev_lower[i], diff_dev_upper[i]) ,
        rep(n + 2 - i - 0.5, 2) ,
        lwd = 0.5 ,
        col = dcol
      )
    }
  }
})


setMethod( "plot" , "precis" , function(x,y,...) precis_plot(x,y,...) )
setMethod( "plot" , "coeftab" , function(x,y,...) coeftab_plot(x,y,...) )
