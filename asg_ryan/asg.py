#Name: Nino Paul Batanay

#Prob 6.13
def logbx_secant(x0, x1, b, x, iters=1000, eps=1e-10):
    def ln_taylor(x, n=10000):
        return 2*sum([1/(2*i+1)*((x-1)/(x+1))**(2*i+1) for i in range(n)])
    f = lambda y,b,x: y - ln_taylor(x)/ln_taylor(b)

    x_0 = x0
    x_1 = x1
    for i in range(iters):
        x_2 = x_1 - f(x_1, b, x)*(x_1-x_0)/(f(x_1, b, x)-f(x_0, b, x))
        if abs(x_2-x_1) < abs(x_1)*eps:
            break
        x_0,x_1 = x_1,x_2
    return x_2

#Prob 6.14
def newton(f, x0, f_prime, n_iters=100, eps=1e-10):
    xold = x0
    xnew = None
    for i in range(n_iters):
        xnew = xold - f(xold)/f_prime(xold)
        if abs(xnew-xold) < abs(xold)*eps:
            break
        xold = xnew
    return xnew

#Prob 6.15
def secant(x0, x1, epsilon, f_prime, iters=1000):
    x2 = -np.inf
    ind = True
    count = 0
    while ind and count<iters:
        x2 = x1 - f_prime(x1)*(x1 - x0)/(f_prime(x1) - f_prime(x0))
        ind = np.abs(x2-x1)>np.abs(x1)*epsilon
        x0,x1 = x1,x2
        count+=1
    return x2
