module noise
using Random


export G,select_gradient,ψ,perlin_mesh

#a set of unit vectors to select gradients from
const G::Vector{Vector{Float64}} = [
    [1,1]./ (√2),
    [-1,1]./ (√2),
    [1,-1]./ (√2),
    [-1,-1]./ (√2),
    [1.,0.],
    [0.,1.],
    [-1.,0.],
    [0.,-1.]
]

"""
Calculates pseudo randomly selected gradient vectors on the edges of a lattice. Lattice is a xn by yn grid tiled by l.

Args:
    xn::Int64 - The width of the mesh
    yn::Int64 - The length of the mesh
    l::Int64 - The size of the lattice to overlay on the mesh
    seed::Int64 - The seed that determines how to shuffle gradient vectors across the mesh

Returns:
    Array{Float64,3} - An [xn,yn,2] array where entry [i,j,:] is the random gradient vector assigned to the lattic at corner i,j

"""
function select_gradient(xn::Int64,yn::Int64,l::Int64,seed::Int64)::Array{Float64, 3}
    #create a permutation table to choose enough gradient vector
    xp::Int64 = (xn/l)+1
    yp::Int64 = (yn/l)+1
    Random.seed!(seed)
    perm = shuffle(0:(xp*yp)-1)
    #now to choose vectors at latice coordinates we just pick out G[i%8]
    grads = zeros(xp,yp,2)
    k = 1
    for i=1:xp,j=1:yp
        grads[i,j,:] = G[ (perm[k]%8)+1 ]
        k+=1
    end
    #now accesing grads[i,j,:]
    grads
end

"""
Maps a value on a range [x0,x1], to a range of [0,1].

Args: 
    x0 - the lower range value
    x1 - the higher range value
    r - the actual value to scale
Returns:
    The value of r when scaled to be in the range of [0,1].
"""
function mapval(x0,x1,r)
    slope = 1/(x1-x0)
     slope*(r - x0)
end

"""
Perlin's improved fading function. Maps the unit square [0,1]x[0,1] to the range [0,1].
"""
ψ(x,y) = (6x^5 - 15x^4 + 10x^3)*(6y^5 - 15y^4 + 10y^3)


"""
Calculate an [xn,yn] matrix that represents perlin noise on a 2D mesh.

Args:
    xn::Int64 - The width of the mesh
    yn::Int64 - The length of the mesh
    l::Int64 - The size of the lattice to overlay on the mesh
    p::Float64 - A value that determines how fast additive noise terms die off
    k::Int64 - How many additional noise terms to sum
    seed::Int64 = 10 - The seed that determines how to shuffle gradient vectors across the mesh 
"""
function perlin_mesh(xn::Int64,yn::Int64,l::Int64,p::Float64,k::Int64;seed::Int64 = 10)::Matrix{Float64}
    #correct xn and yn to make sure the whole grid is l^2 tile-able
    xn = round(Int64,xn/l)*l
    yn = round(Int64,yn/l)*l

    #select a random gradient vector at each lattice point
    g = select_gradient(xn,yn,l,seed)

    #function that maps an index of a lattice to the actual position
    f(x) = (x .- 1) .* l

    #create an array storing
    # 2 entries for a position, 8 entries for the 4 gradients, 8 entries for 4 displacments from gradients
    M = zeros(xn,yn,18)

    for i=1:xn,j=1:yn
        #place the positions
        x = i - 0.5
        y = j - 0.5
        r = [x y]
        M[i,j,1:2] = r
        #calculate the corner indices
        x0::Int64,y0::Int64 = floor( (i-1) /l)+1,floor( (j-1) /l)+1
        x1,y1 = x0+1,y0+1
        #place the 4 gradient vectors
        M[i,j,3:10] = [g[x0,y0,:]... g[x1,y0,:]... g[x0,y1,:]... g[x1,y1,:]... ]
        #get the coordinates of the corners
        x_0,y_0 = f([x0 y0])
        x_1,y_1 = f([x1 y1])
        #construct the weights by mapping the position onto a [0,1] range then creating the fade function weights
        modx = mapval(x_0,x_1,x)
        mody = mapval(y_0,y_1,y)
        fade_pre = [ψ(1-modx,1-mody) ψ(modx,1-mody) ψ(1-modx,mody) ψ(modx,mody)]
        fades = reshape([fade_pre;fade_pre],1,2*size(fade_pre,2))
        #place the displacements from the coordinates weighted by fade values
        M[i,j,11:18] = fades .* [(r-[x_0 y_0])... (r-[x_1 y_0])... (r-[x_0 y_1])... (r-[x_1 y_1])...]
    end

    #multiply the entries of the gradients times the entries of the wieghted displacements
    #then sum them across the 3 dimension to get the noise at each point
    noise = dropdims(sum(M[:,:,3:10] .* M[:,:,11:end],dims=3),dims=3)

    #make a copy of the initial noise turn and start summing octaves oof it
    perlin = copy(noise)
    for i in 1:k
        #shuffle the indices accordingly
        x_shuff = ( ( 2^i .* (0:xn-1) ) .% xn ) .+ 1
        y_shuff = ( ( 2^i .* (0:yn-1) ) .% yn ) .+ 1
        #calculate the new noise term
        perlin += p^i .* noise[x_shuff,y_shuff]
    end

    #rescale each value so it lies between 0 and 1.0
    p_max = maximum(perlin)
    p_min = minimum(perlin)
    perlin_scale(x) = mapval(p_max,p_min,x)
    perlin = perlin_scale.(perlin)
    #return the noise mesh
    return perlin
end


end;