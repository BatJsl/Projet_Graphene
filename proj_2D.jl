using Plots
using FFTW
using KrylovKit

#Nos paramétre
struct parameters_2D
    Kmax
    σ
    Lx
    Ly
end

P_2D = parameters_2D(100.,1.,1.,1.)

function Vpotentiel(x,y)
    -exp(-(x-0.5)^2/P_2D.σ)*exp(-(y-0.5)^2/P_2D.σ)
end

function f(x,y)
    cos(2*pi*x)+sin(2*pi*y)
end

function evaluation(f)
    eval = zeros((Int64(2*P_2D.Kmax)),(Int64(2*P_2D.Kmax)))
    for i in 1:(Int64(2*P_2D.Kmax))
        for j in 1:(Int64(2*P_2D.Kmax))
            eval[i,j] = f((j-1)/(2*P_2D.Kmax),(i-1)/(2*P_2D.Kmax))
        end
    end
    eval
end

V = evaluation(Vpotentiel)
kx =  P_2D.Kmax*2 .* fftfreq(Int64(2*P_2D.Kmax))

Ii = zeros(Int64(2*P_2D.Kmax),Int64(2*P_2D.Kmax))
for i in 1:Int64(2*P_2D.Kmax)
    for j in 1:Int64(2*P_2D.Kmax)
        Ii[i,j]=kx[i]^2+kx[j]^2
    end
end

a = evaluation(f)

Δ = (X -> (4*pi^2 .* X .* Ii))

CFV = (X -> fft(ifft(X) .* V))

H = (X -> Δ(X) + CFV(X))

ii = ones((Int64(2*P_2D.Kmax)),(Int64(2*P_2D.Kmax)))

#La résolution
res = eigsolve(H,ii,4,:SR)

(E,VP) = res

println(E)

x = [i*P_2D.Lx / (2*P_2D.Kmax) for i in 0:(2*Int64(P_2D.Kmax)-1)]
y = [i*P_2D.Lx / (2*P_2D.Kmax) for i in 0:(2*Int64(P_2D.Kmax)-1)]
#plot(0:(2*Int64(P.Kmax) - 1), real.(ifft(cos.(2*pi .*x))) )
surface(x,y,real.(conj.(ifft(VP[1])) .* ifft(VP[1])))
