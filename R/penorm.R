penorm <-
function (q, m = 1, sd = 2) 
{
    x = (q - m)/sd
    p = pnorm(x)
    d = dnorm(x)
    u = -d - z * p
    asy = u/(3 * u + x)
    return(asy)
}
