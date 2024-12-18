### A Pluto.jl notebook ###
# v0.19.40

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ a9208590-a822-11ef-04a1-735c2a12a98a
using PlutoUI

# ╔═╡ 1feb3ab7-a77f-454a-aaa7-9fdb087e680b
begin
	using Plots
	using LinearAlgebra
	using Statistics
end

# ╔═╡ 2322aa25-07bf-46c2-9e74-6ecf2255026a
PlutoUI.TableOfContents(title="Nonlinear FEM: Application", aside=true)

# ╔═╡ 8587302a-08ed-4add-88be-8df3e334c539
md"""
Notebook done by Juan Carlos Galvis and Carlos Nosa.
"""

# ╔═╡ d8f6f94a-a8d3-4f22-a0eb-e4a144f47e37
md"""
We are going to use the following libraries:
"""

# ╔═╡ 6d2ae9c5-ca56-4b31-ab70-6c26eb187a81
md"""
# Theoretical part
"""

# ╔═╡ 6f14c2e3-33e9-4ad0-89d7-ed3ea7e8dc45
md"""
**Strong formulation**

$(S)\begin{cases}
\phi''(t) = A(\psi(t)+\psi_0)\ e^{\gamma \phi(t)},\\
\psi''(t) = -\frac{A}{\delta}\ e^{\gamma \phi(t)},\\
\phi(0) = \phi(2\pi),  \phi'(0)=\phi'(2\pi),\\
\psi(0) = \psi(2\pi),  \psi'(0)=\psi'(2\pi).\\
\end{cases}$
"""

# ╔═╡ e9e31182-647b-437d-ac04-eab584a10c5d
md"""
**Strong formulation for Picard**

$(S)\begin{cases}
-\phi''(t) +\phi(t) = -A(\psi(t)+\psi_0)\ e^{\gamma \phi(t)}+\phi(t),\\
-\psi''(t)+\psi(t) = \frac{A}{\delta}\ e^{\gamma \phi(t)}+\psi,\\
\phi(0) = \phi(2\pi),  \phi'(0)=\phi'(2\pi),\\
\psi(0) = \psi(2\pi),  \psi'(0)=\psi'(2\pi).\\
\end{cases}$
"""

# ╔═╡ 89763b3b-69c4-40da-b060-7eeb0e4b86f7


# ╔═╡ b8bba7a2-c33d-4e7d-b941-f6c4074b7c23
md"""
**Weak formulation for Picard**

$(W)\begin{cases}
\displaystyle\int_{0}^{2\pi}\phi'p'+\phi p = -A\int_{0}^{2\pi} (\psi+\psi_0)\ e^{\gamma \phi}p+\phi p\\
\displaystyle\int_{0}^{2\pi}\psi' q' +\psi q = \frac{A}{\delta}\int_{0}^{2\pi} e^{\gamma \phi}q 
+\psi q\end{cases}$
"""

# ╔═╡ 733e9c4f-d122-46a9-a908-14e2fc4867ea


# ╔═╡ 22099ce3-7655-4652-8e6c-540de4c61548
md"""
# Implementation
"""

# ╔═╡ 85e8541b-b1dc-423b-b8f2-a3dee5e51992
md"""
## Picard's method
"""

# ╔═╡ 05c8fb62-ba89-46c5-b182-8999ea519d06
md"""
Given an initial aproximation $\psi_i$ y $\phi_i$, we build the sequences $\{\psi_{n}\}$ and $\{\phi_{n}\}$ in the following way: Given $\psi_n$ and $\phi_n$, the functions $\psi_{n+1}$ and $\phi_{n+1}$ are defined as the functions satisfying the following equations for all test functions $p$ and $q$:
"""

# ╔═╡ 89414081-c3ad-42b8-a409-166eac20764f
md"""
**Picard's formulation**

$(P)\begin{cases}
\displaystyle\int_{0}^{2\pi}\phi_{n+1}'p'+\phi_{n+1} p = -A\int_{0}^{2\pi} (\psi_{n}+\psi_0)\ e^{\gamma \phi_{n}}p+\phi p\\
\displaystyle\int_{0}^{2\pi}\psi{n+1}' q' +\psi{n+1} q = \frac{A}{\delta}\int_{0}^{2\pi} e^{\gamma \phi_{n}}q 
+\psi_{n} q\end{cases}$
"""

# ╔═╡ 2269f1c8-26e3-4426-8037-756566fb9d63
md"""
If we represent the unknown functions as

$\phi_{n+1}(t) = \sum_{j=0}^{M}h_{j}b_{j}(t)$

$\psi_{n+1}(t) = \sum_{j=0}^{M}k_{j}b_{j}(t)$

where $M$ is a positive integer, $h_0 = h_M$ and $k_0 = k_M$, we can write the previous equations as

$D h_{n+1} =b_{n}$

$D k_{n+1} = c_{n}$

where
"""

# ╔═╡ 53772983-f2dc-4b93-99f6-eada9726c2c5
function SolNum1Picard(M,maxiter,ψ₀,A,δ,γ)
	#Pesos y puntos de integración numérica en [-1,1]
	ω_r = [5/9,8/9, 5/9]
	ζ_r = [-sqrt(15)/5,0,+sqrt(15)/5]

	#Matrices resultantes de FEM
	
	D = zeros(M,M)
	E = zeros(M,M)
	b = zeros(M,1)
	c = zeros(M,1)

	# k = ϕᵢ * ones(M+1,1)
	# h = ψᵢ * ones(M+1,1)

	h=zeros(maxiter,M)
	k=zeros(maxiter,M)
	nodes=0:2*π/M:2*π;
	k[1,:] = sin.(nodes[1:M])
	h[1,:] = cos.(nodes[1:M])

	
	for iter in 2:maxiter
	
	for i=1:M #Iteración por cada elemento
		
		#Definición del elemento
		x0=2*π*(i-1)/M
		x1=2*π*i/M
		Δ=x1-x0
#	f = [exp(γ* ((h[i+1]-h[i])*t/4 + h[i]) )+(h[i+1]-h[i])*t/4 + h[i] for t in 1:3]
		#Pesos y puntos de cuadratura locales
		ω = 0.5*(x1-x0)*ω_r
		ζ = [x0,x0,x0]+ 0.5*([1,1,1]+ζ_r)*(Δ) 
		
		#Funciones base locales
		ϕ1  = ([x1,x1,x1]-ζ)/Δ 
		ϕ2  = (ζ-[x0,x0,x0])/Δ 
		∂ϕ1 = [-1,-1,-1]/Δ
		∂ϕ2 = [1,1,1]/Δ

		#Grados de libertad locales
		dof = [i,i+1] 
		if i == M
			dof = [M,1]
		end

		hζ= ϕ1*h[iter-1,dof[1]]+ϕ2*h[iter-1,dof[2]]
		kζ= ϕ1*k[iter-1,dof[1]]+ϕ2*k[iter-1,dof[2]]
		f=-A*(kζ .+ψ₀).*exp.(γ*hζ) +hζ
		g=(A/δ)*exp.(γ*hζ) .+kζ
		#Submatrices locales
		d11  = sum(ω.*∂ϕ1.*∂ϕ1 + ω.*ϕ1.*ϕ1)
	    d12  = sum(ω.*∂ϕ1.*∂ϕ2 + ω.*ϕ1.*ϕ2)
	    d21  = sum(ω.*∂ϕ2.*∂ϕ1 + ω.*ϕ2.*ϕ1)
	    d22  = sum(ω.*∂ϕ2.*∂ϕ2 + ω.*ϕ2.*ϕ2)
		Dloc = [d11 d12; d21 d22]

#		e11  = sum(ω.*f.*ϕ1.*ϕ1)
#	    e12  = sum(ω.*f.*ϕ1.*ϕ2)
#	    e21  = sum(ω.*f.*ϕ2.*ϕ1)
#	    e22  = sum(ω.*f.*ϕ2.*ϕ2)
#		Eloc = [e11 e12; e21 e22]
		
	    b1   = sum(ω.*f.*ϕ1)
	    b2   = sum(ω.*f.*ϕ2)
		bloc = [b1,b2]
	    c1   = sum(ω.*g.*ϕ1)
	    c2   = sum(ω.*g.*ϕ2)
		cloc = [c1,c2]

		#Ensamblamiento de matrices globales   
		D[dof,dof]  = D[dof,dof] + Dloc
#		E[dof,dof]  = E[dof,dof] + Eloc
		b[dof]      = b[dof] + bloc	 
		c[dof]      = c[dof] + cloc	 
	end
	k[iter,1:M] = D[1:M,1:M] \ (c[1:M])
	h[iter,1:M] = D[1:M,1:M] \ (b[1:M])


	end
	return(h[maxiter,:],k[maxiter,:])
end

# ╔═╡ 097d2cec-1f08-4f45-9112-04a1687230ed
begin
	ψ₀ = 1
	A = 1
	δ = 1000
	γmin=0
	γmax=5
end

# ╔═╡ 7d503c27-d4e2-4ea5-b960-4690a03932be
begin
	iter = 100
	M = 100
	nodes=0:2*π/M:2*π
	nodesp=nodes[1:M]

	γrange=γmin:.1:γmax;
	solhγ=zeros(length(γrange),M)
	solkγ=zeros(length(γrange),M)
	for i=1:length(γrange)
	hsol,ksol=SolNum1Picard(M,iter,ψ₀,A,δ,γrange[i] )
	solhγ[i,:]=hsol
	solkγ[i,:]=ksol
	end
	
end

# ╔═╡ 865cb63f-3d65-4226-ac9c-6204e64829e4
	@bind iγ Slider(1:length(γrange), default=1)

# ╔═╡ a7fd9bcc-33d0-4f3b-8022-6644eb1127a2
	plot(plot(nodesp,solhγ[iγ,:],label="ϕ", lw=4, ylims=(-10, 10)),
		plot(nodesp,solkγ[iγ,:],label="ψ", lw=4,ylims=(-10, 10)),
		layout=(1,2))

# ╔═╡ 2c1d3476-ab68-4a17-b783-9b5c1f200558
# ╠═╡ disabled = true
#=╠═╡
anim = @animate  for i in 1:iter
	plot(plot(nodesp,hsol[i,:],label="ϕ", lw=4, ylims=(-1.5, -0.5)),
		plot(nodesp,ksol[i,:],label="ψ", lw=4,ylims=(0.5, 1)),
		layout=(1,2), title="Picard(Iteration: $i / $iter)",)
end
  ╠═╡ =#

# ╔═╡ d6e4b616-157b-4694-8385-55127fb482be
#=╠═╡
# Save the animation as a GIF
gif(anim, "Picard2.gif", fps=15)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─a9208590-a822-11ef-04a1-735c2a12a98a
# ╟─2322aa25-07bf-46c2-9e74-6ecf2255026a
# ╟─8587302a-08ed-4add-88be-8df3e334c539
# ╟─d8f6f94a-a8d3-4f22-a0eb-e4a144f47e37
# ╠═1feb3ab7-a77f-454a-aaa7-9fdb087e680b
# ╟─6d2ae9c5-ca56-4b31-ab70-6c26eb187a81
# ╟─6f14c2e3-33e9-4ad0-89d7-ed3ea7e8dc45
# ╟─e9e31182-647b-437d-ac04-eab584a10c5d
# ╠═89763b3b-69c4-40da-b060-7eeb0e4b86f7
# ╠═b8bba7a2-c33d-4e7d-b941-f6c4074b7c23
# ╠═733e9c4f-d122-46a9-a908-14e2fc4867ea
# ╟─22099ce3-7655-4652-8e6c-540de4c61548
# ╟─85e8541b-b1dc-423b-b8f2-a3dee5e51992
# ╟─05c8fb62-ba89-46c5-b182-8999ea519d06
# ╠═89414081-c3ad-42b8-a409-166eac20764f
# ╠═2269f1c8-26e3-4426-8037-756566fb9d63
# ╠═53772983-f2dc-4b93-99f6-eada9726c2c5
# ╠═097d2cec-1f08-4f45-9112-04a1687230ed
# ╠═7d503c27-d4e2-4ea5-b960-4690a03932be
# ╠═865cb63f-3d65-4226-ac9c-6204e64829e4
# ╠═a7fd9bcc-33d0-4f3b-8022-6644eb1127a2
# ╠═2c1d3476-ab68-4a17-b783-9b5c1f200558
# ╠═d6e4b616-157b-4694-8385-55127fb482be
