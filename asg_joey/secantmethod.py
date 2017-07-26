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
