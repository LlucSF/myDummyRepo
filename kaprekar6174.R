#################EXTRA FOOTAGE#########################


numb <-0:999
max1digit <- 0:999
max2digit <- 0:999
max3digit <- 0:999
max4digit <- 0:999

for (i in 1:1000) {
  max4digit[i] <- sort(extractDigits(numb[i]))[1]
}

for (i in 1:1000) {
  max3digit[i] <- sort(extractDigits(numb[i]))[2]
}

for (i in 1:1000) {
  max2digit[i] <- sort(extractDigits(numb[i]))[3]
}

for (i in 1:1000) {
  max1digit[i] <- sort(extractDigits(numb[i]))[4]
}

data <- data.frame(digit = c(max1digit,max2digit,max3digit,max4digit),
                   func = rep(c("First","Second","Third","Fourth"),each = 1000),
                   numb = rep(numb,times = 4))


library(ggplot2)
ggplot() + geom_path(data =data,mapping = aes(x = numb, y = digit, color = func), size = 1.2) + theme_bw()












### setting variables and functions ###
numbers <- 0:9999


extractDigits <- function(number)
{
  digits <- c(trunc(number/1000),
              (trunc(number/100)-trunc(number/1000)*10),
              (trunc(number/10)-(trunc(number/100)*10)),
              trunc(number)-(trunc(number/1000)*1000)-100*(trunc(number/100)-trunc(number/1000)*10)-10*(trunc(number/10)-(trunc(number/100)*10))
              )
  
return(digits)
  
}

pasteDigits <- function(digits)
{
  return(digits[1]*1000+digits[2]*100+digits[3]*10+digits[4])
}



### keprak loop ### 
{
  iterations <- rep(0,times = 10000)
  iterations[6175] <- 1
  start_time <- Sys.time()
  for(i in numbers)
  {
    digits <- extractDigits(i)
    if(!(length(unique(digits))==1)) ##Avoid number 6174 and same digits numbers (0000,1111,...)
    {

      while(!all(digits == c(6,1,7,4))) ##while digits are nor 6 1 7 4
      {
        digits <- extractDigits(pasteDigits(sort(digits,decreasing = T))-pasteDigits(sort(digits)))
        iterations[i+1] <- iterations[i+1] + 1
        if(pasteDigits(digits) < i)
        {
          iterations[i+1] <-  iterations[i+1] + iterations[pasteDigits(digits)+1]
          break
        }
      }
    }
  }
  end_time <- Sys.time()
  print(end_time-start_time)
}



### plot ###
ggplot2::ggplot() +
  ggplot2::geom_point(mapping = ggplot2::aes(x = numbers, y = iterations),size = 0.6) + 
  ggplot2::theme_bw()

ggplot2::ggplot() +
  ggplot2::geom_histogram(mapping = ggplot2::aes( x = iterations),size = 0.6,bins = 8) + 
  ggplot2::theme_bw()


iteraionImg <- matrix(iterations, nrow = 100, ncol = 100)
par(mar = c(0, 0, 0, 0), mgp = c(0, 0, 0), mfrow = c(1,1), las = 0, mai = c(0,0,0,0))

img <- raster::raster(nrows = nrow(iteraionImg),
                      ncols = ncol(iteraionImg),
                      xmn = 0,
                      xmx = ncol(iteraionImg),
                      ymn = 0,
                      ymx = nrow(iteraionImg)
) 

raster::values(img) <- iteraionImg 

raster::plot(x = img,
             col =rainbow(7),
             useRaster = T,
             interpolate = T,
             box = F,
             bty = "n",
             axes = F,
             legend = F,
             
)


## 3D plot ##  ##5063 ##2538

iterations[which(iterations==0)] <- 14
iteraionImg <- matrix(iterations, nrow = 100, ncol = 100)

plotly::plot_ly(z = iteraionImg, colors = rainbow(8)) %>% plotly::add_surface(showlegend = F)  


smth <- fields::interp.surface.grid(obj = list(x = 1:100, y = 1:100, z = iteraionImg),
                                    grid.list = list(x = seq(1,100,length.out = 1000),y =seq(1,100,length.out = 1000))
                                    )

plotly::plot_ly(z = ~smth$z) %>% plotly::add_surface(showlegend = F)  






