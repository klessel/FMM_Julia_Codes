function svd_ranktol(A::Array{Float64,2},r::Int64)
#INPUT:             A       Array{Float64}(2)
#                   r       (Uint64) largest allowable rank of hankel blocks -
#                           this determines the amount of compression, and
#                           should be compatible with your input tree.
#                           (if the maximum allowable rank is p, then the
#                           tree partitions should be no smaller than 3p)
#                           Corresponding singular vectors below this rank
#                           will be dropped.
#OUTPUT:            Uhat    (Array{Float64}(2)) left singular vectors of A
#                           corresponding to singular values higher than
#                           the chosen tolerance
#                   Vhat    (Array{Float64}(2)) right singular vectors of A
#                           corresponding to singular values higher than
#                           the chosen tolerance
#                   Shat    (Array{Float64}(1)) array containing the singular
#                           values of A higher than the chosen tolerance
    U, S, V = svd(A);
    if isempty(A)
        # # julia handles svd's of empty matrices incorrectly at this time.
        # These if statements correct that
        m, n = size(A);
        if size(A,2) == 0
            Uhat = Matrix{Float64}(I,m,n); #eye(m,n);
            Shat = Array{Float64,1}(undef,0); #for some reason the function diagm() (used in the fmm construction code) will only accept this input to produce a 0x0 matrix, and errors otherwise
            Vhat = Array{Float64,2}(undef,0,0);
        elseif size(A,1) == 0
            Uhat = Array{Float64,2}(undef,0,0);
            Shat = Array{Float64,1}(undef,0);
            Vhat = Matrix{Float64}(I,n,m);#eye(n,m);
        end
    else
        Uhat = U[:,1:r];
        Shat = S[1:r];
        Vhat = V[:,1:r];
    end
    Uhat, Shat, Vhat
 end
