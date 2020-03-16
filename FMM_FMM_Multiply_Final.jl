#
# Author: Kristen Lessel, kristenlessel@gmail.com
# June 15 2018
# FMM x FMM Multiply and utilities
# This version allows for passage of Arrays of Arrays or OffsetArrays of Arrays
# This version stores F, B anc B_C as Offset Arrays instead of using functions

function FMM_FMM_Multiply(D::OffsetArray{OffsetArray{Array{Float64}}},U::OffsetArray{Array{Float64}},R::Array{OffsetArray{Array{Float64}}},B::Array{OffsetArray{OffsetArray{Array{Float64}}}},W::Array{OffsetArray{Array{Float64}}},V::OffsetArray{Array{Float64}},D̃::OffsetArray{OffsetArray{Array{Float64}}},Ũ::OffsetArray{Array{Float64}},R̃::Array{OffsetArray{Array{Float64}}},B̃::Array{OffsetArray{OffsetArray{Array{Float64}}}},W̃::Array{OffsetArray{Array{Float64}}},Ṽ::OffsetArray{Array{Float64}},depth::Int64)
    #Multiplies 2 fmm matrices, with U,R,B,W,V being from the first 3 point fmm structure,
    #and Ũ,R̃,B̃,W̃,Ṽ being from the second 3 point fmm structure. Output is a 5 point fmm
    #structure stored in arrays as listed here:
    #Input:             D:: OffsetArray{OffsetArray{Array{Float64}}}
    #                   U,V:: OffsetArray{Array{Float64}}
    #                   R,W:: Array{OffsetArray{Array{Float64}}}
    #                   B:: Array{OffsetArray{Array{Array{Float64}}}}
    #                   depth:: Int64  (depth of the two input fmm structures)
    #Output:            D_C:: OffsetArray{OffsetArray{Array{Float64}}}
    #                   U_C,V_C:: Array{Array{Float64}}
    #                   R_C,W_C:: Array{Array{Array{Float64}}}
    #                   B_C:: Array{Array{Array{Float64}}}

    ################################################################################
    #Compute U_C, V_C
    U_C = Array{Array{Float64}}(undef,2^(depth)+1);
    V_C = Array{Array{Float64}}(undef,2^(depth)+1);
    for i = 1:2 .^depth
        U_C[i] = [U[i] D[i][i-1]*Ũ[i-1] D[i][i]*Ũ[i] D[i][i+1]*Ũ[i+1]];
        V_C[i] = [Ṽ[i] D̃[i-1][i]'*V[i-1] D̃[i][i]'*V[i] D̃[i+1][i]'V[i+1]];
    end

    ################################################################################
    # Upsweep to Compute G
    G = Array{OffsetArray{Array{Float64}}}(undef,depth);
    for k = 1:depth
        G[k] = OffsetArray{Array{Float64}}(undef,-1:2^(k)+4);
    end

    #Compute G at the leaves
    for i = 1:2^depth
        G[depth][i] = V[i]'*Ũ[i];
    end

    #Compute G at the nodes
    for k = depth:-1:2
        for i = 1:2^(k-1)
            G[k-1][i] = W[k][2i-1]'*G[k][2i-1]*R̃[k][2i-1]+W[k][2i]'*G[k][2i]R̃[k][2i];
        end
    end

    #assign empty G matrices to the right and left of the edge of each tree
    #(matrices outside the indices 1:2^k are 0x0 empty matrices)
    for k = 1:depth
        G[k][-1] = zeros(0,0);
        G[k][0] = zeros(0,0);
        G[k][2^k+1] = zeros(0,0);
        G[k][2^k+2] = zeros(0,0);
        G[k][2^k+3] = zeros(0,0);
        G[k][2^k+4] = zeros(0,0);
    end

    ################################################################################
    #Compute R_C, W_C
    R_C = Array{Array{Array{Float64}}}(undef,depth);
    W_C = Array{Array{Array{Float64}}}(undef,depth);
    for k = 1:depth
        R_C[k] = Array{Array{Float64}}(undef,2^(k)+2);
        W_C[k] = Array{Array{Float64}}(undef,2^(k)+2);
    end

    R_C[1][1] = R[1][1]; #if i wanted to add this into the loop I would just have to store the 4 empty columns of matrices of B1;*,*
    R_C[1][2] = R[1][2];
    W_C[1][1] = W[1][1];
    W_C[1][2] = W[1][2];

    for k = 2:depth

        for i = 1:2^(k-1)

            #create zero blocks of the appropriate dimensions
            ml_Rl, nl_Rl = size(R[k][2i-3]); #m_Rl = m_Rr
            ml_Rr, nl_Rr = size(R[k][2i-2]);
            m0_Rl, n0_Rl = size(R[k][2i-1]);
            m0_Rr, n0_Rr = size(R[k][2i]);
            mr_Rl, nr_Rl = size(R[k][2i+1]);
            mr_Rr, nr_Rr = size(R[k][2i+2]);

            ml_R̃l, nl_R̃l = size(R̃[k][2i-3]);
            ml_R̃r, nl_R̃r = size(R̃[k][2i-2]);
            m0_R̃l, n0_R̃l = size(R̃[k][2i-1]);
            m0_R̃r, n0_R̃r = size(R̃[k][2i]);
            mr_R̃l, nr_R̃l = size(R̃[k][2i+1]);
            mr_R̃r, nr_R̃r = size(R̃[k][2i+2]);

            ml_Wl, nl_Wl = size(W[k][2i-3]); #m_Wl = m_Wr
            ml_Wr, nl_Wr = size(W[k][2i-2]);
            m0_Wl, n0_Wl = size(W[k][2i-1]);
            m0_Wr, n0_Wr = size(W[k][2i]);
            mr_Wl, nr_Wl = size(W[k][2i+1]);
            mr_Wr, nr_Wr = size(W[k][2i+2]);

            ml_W̃l, nl_W̃l = size(W̃[k][2i-3]);
            ml_W̃r, nl_W̃r = size(W̃[k][2i-2]);
            m0_W̃l, n0_W̃l = size(W̃[k][2i-1]);
            m0_W̃r, n0_W̃r = size(W̃[k][2i]);
            mr_W̃l, nr_W̃l = size(W̃[k][2i+1]);
            mr_W̃r, nr_W̃r = size(W̃[k][2i+2]);

            mrr_R̃l, nrr_R̃l =size(R̃[k][2i+3]);
            mrr_Wl, nrr_Wl =size(W[k][2i+3]);

            R_C[k][2i-1] = [R[k][2i-1] B[k][2i-1][2i-3]*G[k][2i-3]*R̃[k][2i-3] zeros(m0_Rl,n0_R̃l) B[k][2i-1][2i+1]*G[k][2i+1]*R̃[k][2i+1]+B[k][2i-1][2i+2]*G[k][2i+2]*R̃[k][2i+2];
                            zeros(ml_R̃r,n0_Rl) R̃[k][2i-2] zeros(ml_R̃r,n0_R̃l) zeros(ml_R̃r,nr_R̃l);
                            zeros(m0_R̃l,n0_Rl) zeros(m0_R̃l,nl_R̃r) R̃[k][2i-1] zeros(m0_R̃l,nr_R̃l);
                            zeros(m0_R̃r,n0_Rl) zeros(m0_R̃r,nl_R̃r) R̃[k][2i] zeros(m0_R̃r,nr_R̃l)];
            W_C[k][2i-1] = [W̃[k][2i-1] B̃[k][2i-3][2i-1]'*G[k][2i-3]'*W[k][2i-3] zeros(m0_W̃l,n0_Wl) B̃[k][2i+1][2i-1]'*G[k][2i+1]'*W[k][2i+1] + B̃[k][2i+2][2i-1]'*G[k][2i+2]'*W[k][2i+2];
                            zeros(ml_Wr,n0_W̃l) W[k][2i-2] zeros(ml_Wr,n0_Wl) zeros(ml_Wr,nr_Wl);
                            zeros(m0_Wl,n0_W̃l) zeros(m0_Wl,nl_Wr) W[k][2i-1] zeros(m0_Wl,nr_Wl);
                            zeros(m0_Wr,n0_W̃l) zeros(m0_Wr,nl_Wr) W[k][2i] zeros(m0_Wr,nr_Wl)];

            R_C[k][2i] = [R[k][2i] B[k][2i][2i-3]*G[k][2i-3]*R̃[k][2i-3] + B[k][2i][2i-2]*G[k][2i-2]*R̃[k][2i-2] zeros(m0_Rr,n0_R̃l) B[k][2i][2i+2]*G[k][2i+2]*R̃[k][2i+2];
                          zeros(m0_R̃l,n0_Rr) zeros(m0_R̃l,nl_R̃l) R̃[k][2i-1] zeros(m0_R̃l,nr_Rr);
                          zeros(m0_R̃r,n0_Rr) zeros(m0_R̃r,nl_R̃l) R̃[k][2i] zeros(m0_R̃r,nr_Rr);
                          zeros(mr_R̃l,n0_Rr) zeros(mr_R̃l,nl_R̃l) zeros(mr_R̃l,n0_R̃l) R̃[k][2i+1]];
            W_C[k][2i] = [W̃[k][2i] B̃[k][2i-3][2i]'*G[k][2i-3]'*W[k][2i-3]+B̃[k][2i-2][2i]'*G[k][2i-2]'*W[k][2i-2] zeros(m0_W̃r,n0_Wl) B̃[k][2i+2][2i]'*G[k][2i+2]'*W[k][2i+2];
                          zeros(m0_Wl,n0_W̃r) zeros(m0_Wl,nl_Wl) W[k][2i-1] zeros(m0_Wl,nr_Wl);
                          zeros(m0_Wr,n0_W̃r) zeros(m0_Wr,nl_Wl) W[k][2i] zeros(m0_Wr,nr_Wl);
                          zeros(mr_Wl,n0_W̃r) zeros(mr_Wl,nl_Wl) zeros(mr_Wl,n0_Wl) W[k][2i+1]];

            if i == 2^(k-1)
                R_C[k][2i+1] = [R[k][2i+1] B[k][2i+1][2i-1]*G[k][2i-1]*R̃[k][2i-1] zeros(mr_Rr,nr_R̃l) B[k][2i+1][2i+3]*G[k][2i+3]*R̃[k][2i+3] + B[k][2i+1][2i+4]*G[k][2i+4]*R̃[k][2i+4];
                                zeros(m0_R̃r,nr_Rl) R̃[k][2i] zeros(m0_R̃r,nr_R̃l) zeros(m0_R̃r,nrr_R̃l);
                                zeros(mr_R̃l,nr_Rl) zeros(mr_R̃l,n0_R̃r) R̃[k][2i+1] zeros(mr_R̃l,nrr_R̃l);
                                zeros(mr_R̃r,nr_Rl) zeros(mr_R̃r,n0_R̃r) R̃[k][2i+2] zeros(mr_R̃r,nrr_R̃l)];

                R_C[k][2i+2] = [R[k][2i+2] B[k][2i+2][2i-1]*G[k][2i-1]*R̃[k][2i-1]+ B[k][2i+2][2i]*G[k][2i]*R̃[k][2i] zeros(mr_Rr,nr_R̃l) B[k][2i+2][2i+4]*G[k][2i+4]*R̃[k][2i+4];
                                zeros(mr_R̃l,nr_Rr) zeros(mr_R̃l,n0_R̃l) R̃[k][2i+1] zeros(mr_R̃l,nrr_R̃l);
                                zeros(mr_R̃r,nr_Rr) zeros(mr_R̃r,n0_R̃l)  R̃[k][2i+2] zeros(mr_R̃r,nrr_R̃l);
                                zeros(mrr_R̃l,nr_Rr) zeros(mrr_R̃l,n0_R̃l) zeros(mrr_R̃l,nr_R̃l)  R̃[k][2i+3]];
                W_C[k][2i+1] = [W̃[k][2i+1] B̃[k][2i-1][2i+1]'*G[k][2i-1]'*W[k][2i-1] zeros(mr_W̃l,nr_Wl) B̃[k][2i+3][2i+1]'*G[k][2i+3]'*W[k][2i+3] + B̃[k][2i+4][2i+1]'*G[k][2i+4]'*W[k][2i+4];
                                zeros(m0_Wr,nr_W̃l) W[k][2i] zeros(m0_Wr,nr_Wl) zeros(m0_Wr,nrr_Wl);
                                zeros(mr_Wl,nr_W̃l) zeros(mr_Wl,n0_Wr) W[k][2i+1] zeros(mr_Wl,nrr_Wl);
                                zeros(mr_Wr,nr_W̃l) zeros(mr_Wr,n0_Wr) W[k][2i+2] zeros(mr_Wr,nrr_Wl)];
                W_C[k][2i+2] = [W̃[k][2i+2] B̃[k][2i-1][2i+2]'*G[k][2i-1]'*W[k][2i-1] zeros(mr_W̃r,nr_Wl) B̃[k][2i+4][2i+2]'*G[k][2i+4]'*W[k][2i+4];
                                zeros(mr_Wl,nr_Wr) zeros(mr_Wl,n0_Wl) W[k][2i+1] zeros(mr_Wl,nrr_Wl);
                                zeros(mr_Wr,nr_Wr) zeros(mr_Wr,n0_Wl) W[k][2i+2] zeros(mr_Wr,nrr_Wl);
                                zeros(mrr_Wl,nr_Wr) zeros(mrr_Wl,n0_Wl) zeros(mrr_Wl,nr_Wl) W[k][2i+3]];
            end

end

end

################################################################################
#Down-sweep recursions (F matrices)

F = Array{OffsetArray{OffsetArray{Array{Float64}}}}(undef,depth);
for k = 1:depth
    F[k] = OffsetArray{OffsetArray{Array{Float64}}}(undef,1:2^k+2);
    for i = 1:2^k+2
        F[k][i] = OffsetArray{Array{Float64}}(undef,i-2:i+2);
    end
end

F[1][1][1] = zeros(0,0);
F[1][1][2] = zeros(0,0);
F[1][1][3] = zeros(0,0);
F[1][2][1] = zeros(0,0);
F[1][3][1] = zeros(0,0);
F[1][2][2] = zeros(0,0);
F[1][2][3] = zeros(0,0);
F[1][3][2] = zeros(0,0);

for k = 2:depth
    for i = 1:2^(k-1)

        #create zero blocks of the appropriate dimensions
        ml_Rl, nl_Rl = size(R[k][2i-3]); #n_Rl = n_Rr
        ml_Rr, nl_Rr = size(R[k][2i-2]);
        m0_Rl, n0_Rl = size(R[k][2i-1]);
        m0_Rr, n0_Rr = size(R[k][2i]);
        mr_Rl, nr_Rl = size(R[k][2i+1]);
        mr_Rr, nr_Rr = size(R[k][2i+2]);

        ml_R̃l, nl_R̃l = size(R̃[k][2i-3]);
        ml_R̃r, nl_R̃r = size(R̃[k][2i-2]);
        m0_R̃l, n0_R̃l = size(R̃[k][2i-1]);
        m0_R̃r, n0_R̃r = size(R̃[k][2i]);
        mr_R̃l, nr_R̃l = size(R̃[k][2i+1]);
        mr_R̃r, nr_R̃r = size(R̃[k][2i+2]);

        ml_Wl, nl_Wl = size(W[k][2i-3]); #m_Wl = m_Wr
        ml_Wr, nl_Wr = size(W[k][2i-2]);
        m0_Wl, n0_Wl = size(W[k][2i-1]);
        m0_Wr, n0_Wr = size(W[k][2i]);
        mr_Wl, nr_Wl =size(W[k][2i+1]);
        mr_Wr, nr_Wr =size(W[k][2i+2]);

        ml_W̃l, nl_W̃l = size(W̃[k][2i-3]);
        ml_W̃r, nl_W̃r = size(W̃[k][2i-2]);
        m0_W̃l, n0_W̃l = size(W̃[k][2i-1]);
        m0_W̃r, n0_W̃r = size(W̃[k][2i]);
        mr_W̃l, nr_W̃l = size(W̃[k][2i+1]);
        mr_W̃r, nr_W̃r = size(W̃[k][2i+2]);

        mrr_R̃l, nrr_R̃l =size(R̃[k][2i+3]);
        mrr_Wl, nrr_Wl =size(W[k][2i+3]);

        #The size of the zero blocks in F correspond to the number of rows in each
        #block of R_C, and the number of columns in each block of W_C^H.

        #F[k][2i-1][2i-1]
        F[k][2i-1][2i-1] = [B[k][2i-1][2i-3]*G[k][2i-3]*B̃[k][2i-3][2i-1] + B[k][2i-1][2i+1]*G[k][2i+1]*B̃[k][2i+1][2i-1] + B[k][2i-1][2i+2]*G[k][2i+2]*B̃[k][2i+2][2i-1] zeros(m0_Rl,ml_Wr) zeros(m0_Rl,m0_Wl) zeros(m0_Rl,m0_Wr);
                            zeros(ml_R̃r,m0_W̃l) zeros(ml_R̃r,ml_Wr) zeros(ml_R̃r,m0_Wl) zeros(ml_R̃r,m0_Wr);
                            zeros(m0_R̃l,m0_W̃l) zeros(m0_R̃l,ml_Wr) zeros(m0_R̃l,m0_Wl) zeros(m0_R̃l,m0_Wr);
                            zeros(m0_R̃r,m0_W̃l) zeros(m0_R̃r,ml_Wr) zeros(m0_R̃r,m0_Wl) zeros(m0_R̃r,m0_Wr)] +
                                R_C[k][2i-1]*F[k-1][i][i]*W_C[k][2i-1]';

        #F[k][2i,2i]
        F[k][2i][2i] = [B[k][2i][2i-3]*G[k][2i-3]*B̃[k][2i-3][2i] + B[k][2i][2i-2]*G[k][2i-2]*B̃[k][2i-2][2i] + B[k][2i][2i+2]*G[k][2i+2]*B̃[k][2i+2][2i] zeros(m0_Rr,m0_Wl) zeros(m0_Rr,m0_Wr) zeros(m0_Rr,mr_Wl);
                        zeros(m0_R̃l,m0_W̃r) zeros(m0_R̃l,m0_Wl) zeros(m0_R̃l,m0_Wr) zeros(m0_R̃l,mr_Wl);
                        zeros(m0_R̃r,m0_W̃r) zeros(m0_R̃r,m0_Wl) zeros(m0_R̃r,m0_Wr) zeros(m0_R̃r,mr_Wl);
                        zeros(mr_R̃l,m0_W̃r) zeros(mr_R̃l,m0_Wl) zeros(mr_R̃l,m0_Wr) zeros(mr_R̃l,mr_Wl)] +
                            R_C[k][2i]*F[k-1][i][i]*W_C[k][2i]';

        #F[k][2i-1,2i]
        F[k][2i-1][2i] = [B[k][2i-1][2i-3]*G[k][2i-3]*B̃[k][2i-3][2i] + B[k][2i-1][2i+2]*G[k][2i+2]*B̃[k][2i+2][2i] zeros(m0_Rl,m0_Wl) zeros(m0_Rl,m0_Wr) B[k][2i-1][2i+1];
                          B̃[k][2i-2][2i] zeros(ml_R̃r,m0_Wl) zeros(ml_R̃r,m0_Wr) zeros(ml_R̃r,mr_Wl);
                          zeros(m0_R̃l,m0_W̃r) zeros(m0_R̃l,m0_Wl) zeros(m0_R̃l,m0_Wr) zeros(m0_R̃l,mr_Wl);
                          zeros(m0_R̃r,m0_W̃r) zeros(m0_R̃r,m0_Wl) zeros(m0_R̃r,m0_Wr) zeros(m0_R̃r,mr_Wl)] +
                              R_C[k][2i-1]*F[k-1][i][i]*W_C[k][2i]';

        #F[k][2i,2i-1]
        F[k][2i][2i-1] = [B[k][2i][2i-3]*G[k][2i-3]*B̃[k][2i-3][2i-1] + B[k][2i][2i+2]*G[k][2i+2]*B̃[k][2i+2][2i-1] B[k][2i][2i-2] zeros(m0_Rr,m0_Wl) zeros(m0_Rr,m0_Wr);
                          zeros(m0_R̃l,m0_W̃l) zeros(m0_R̃l,ml_Wr) zeros(m0_R̃l,m0_Wl) zeros(m0_R̃l,m0_Wr);
                          zeros(m0_R̃r,m0_W̃l) zeros(m0_R̃r,ml_Wr) zeros(m0_R̃r,m0_Wl) zeros(m0_R̃r,m0_Wr);
                          B̃[k][2i+1][2i-1] zeros(mr_R̃l,ml_Wr) zeros(mr_R̃l,m0_Wl) zeros(mr_R̃l,m0_Wr)] +
                              R_C[k][2i]*F[k-1][i][i]*W_C[k][2i-1]';

        #F[k][2i-1,2i+1]
        F[k][2i-1][2i+1] = [zeros(m0_Rl,mr_W̃l) zeros(m0_Rl,m0_Wr) B[k][2i-1][2i+1] B[k][2i-1][2i+2];
                            zeros(ml_R̃r,mr_W̃l) zeros(ml_R̃r,m0_Wr) zeros(ml_R̃r,mr_Wl) zeros(ml_R̃r,mr_Wr);
                            B̃[k][2i-1][2i+1] zeros(m0_R̃l,m0_Wr) zeros(m0_R̃l,mr_Wl) zeros(m0_R̃l,mr_Wr);
                            zeros(m0_R̃r,mr_W̃l) zeros(m0_R̃r,m0_Wr) zeros(m0_R̃r,mr_Wl) zeros(m0_R̃r,mr_Wr)] +
                                R_C[k][2i-1]*F[k-1][i][i+1]*W_C[k][2i+1]';

        #F[k][2i,2i+2]
        F[k][2i][2i+2] = [zeros(m0_Rr,mr_W̃r) zeros(m0_Rr,mr_Wl) B[k][2i][2i+2] zeros(m0_Rr,mrr_Wl);
                          B̃[k][2i-1][2i+2] zeros(m0_R̃l,mr_Wl) zeros(m0_R̃l,mr_Wr) zeros(m0_R̃l,mrr_Wl);
                          B̃[k][2i][2i+2] zeros(m0_R̃r,mr_Wl) zeros(m0_R̃r,mr_Wr) zeros(m0_R̃r,mrr_Wl);
                          zeros(mr_R̃l,mr_W̃r) zeros(mr_R̃l,mr_Wl) zeros(mr_R̃l,mr_Wr) zeros(mr_R̃l,mrr_Wl)] +
                              R_C[k][2i]*F[k-1][i][i+1]*W_C[k][2i+2]';

        #F[k][2i,2i+1]
        F[k][2i][2i+1] = [zeros(m0_Rr,mr_W̃l) zeros(m0_Rr,m0_Wr) zeros(m0_Rr,mr_Wl) B[k][2i][2i+2];
                          B̃[k][2i-1][2i+1] zeros(m0_R̃l,m0_Wr) zeros(m0_R̃l,mr_Wl) zeros(m0_R̃l,mr_Wr);
                          zeros(m0_R̃r,mr_W̃l) zeros(m0_R̃r,m0_Wr) zeros(m0_R̃r,mr_Wl) zeros(m0_R̃r,mr_Wr);
                          zeros(mr_R̃l,mr_W̃l) zeros(mr_R̃l,m0_Wr) zeros(mr_R̃l,mr_Wl) zeros(mr_R̃l,mr_Wr)] +
R_C[k][2i]*F[k-1][i][i+1]*W_C[k][2i+1]';

#F[k][2i+1,2i-1]
F[k][2i+1][2i-1] = [zeros(mr_Rl,m0_W̃l) zeros(mr_Rl,ml_Wr) B[k][2i+1][2i-1] zeros(mr_Rl,m0_Wr);
                    zeros(m0_R̃r,m0_W̃l) zeros(m0_R̃r,ml_Wr) zeros(m0_R̃r,m0_Wl) zeros(m0_R̃r,m0_Wr);
                    B̃[k][2i+1][2i-1] zeros(mr_R̃l,ml_Wr) zeros(mr_R̃l,m0_Wl) zeros(mr_R̃l,m0_Wr);
                    B̃[k][2i+2][2i-1] zeros(mr_R̃r,ml_Wr) zeros(mr_R̃r,m0_Wl) zeros(mr_R̃r,m0_Wr)] +
                        R_C[k][2i+1]*F[k-1][i+1][i]*W_C[k][2i-1]';

#F[k][2i+1,2i]
F[k][2i+1][2i] = [zeros(mr_Rl,m0_W̃r) B[k][2i+1][2i-1] zeros(mr_Rl,m0_Wr) zeros(mr_Rl,mr_Wl);
                  zeros(m0_R̃r,m0_W̃r) zeros(m0_R̃r,m0_Wl) zeros(m0_R̃r,m0_Wr) zeros(m0_R̃r,mr_Wl);
                  zeros(mr_R̃l,m0_W̃r) zeros(mr_R̃l,m0_Wl) zeros(mr_R̃l,m0_Wr) zeros(mr_R̃l,mr_Wl);
                  B̃[k][2i+2][2i] zeros(mr_R̃r,m0_Wl) zeros(mr_R̃r,m0_Wr) zeros(mr_R̃r,mr_Wl)] +
                      R_C[k][2i+1]*F[k-1][i+1][i]*W_C[k][2i]';


#F[k][2i+2,2i]
F[k][2i+2][2i] = [zeros(mr_Rr,m0_W̃r) B[k][2i+2][2i-1] B[k][2i+2][2i] zeros(mr_Rr,mr_Wl);
                  zeros(mr_R̃l,m0_W̃r) zeros(mr_R̃l,m0_Wl) zeros(mr_R̃l,m0_Wr) zeros(mr_R̃l,mr_Wl);
                  B̃[k][2i+2][2i] zeros(mr_R̃r,m0_Wl) zeros(mr_R̃r,m0_Wr) zeros(mr_R̃r,mr_Wl);
                  zeros(mrr_R̃l,m0_W̃r) zeros(mrr_R̃l,m0_Wl) zeros(mrr_R̃l,m0_Wr) zeros(mrr_R̃l,mr_Wl)] +
                      R_C[k][2i+2]*F[k-1][i+1][i]*W_C[k][2i]';


end
end

################################################################################
#Calculation of B's (Expansion Coefficients)

B_C = Array{OffsetArray{OffsetArray{Array{Float64}}}}(undef,depth);
for k = 1:depth
    B_C[k] = OffsetArray{OffsetArray{Array{Float64}}}(undef,1:2^k+5);
    for i = 1:2^k+5
        B_C[k][i] = OffsetArray{Array{Float64}}(undef,i-5:i+5);
    end
end

B_C[1][1][4] = zeros(0,0);
B_C[1][1][5] = zeros(0,0);
B_C[1][1][6] = zeros(0,0);
B_C[1][2][5] = zeros(0,0);
B_C[1][2][6] = zeros(0,0);
B_C[1][4][1] = zeros(0,0);
B_C[1][5][1] = zeros(0,0);
B_C[1][6][1] = zeros(0,0);
B_C[1][5][2] = zeros(0,0);
B_C[1][6][2] = zeros(0,0);

for k = 2:depth
    for i = 1:2^(k-1)-1

        #create zero blocks of the appropriate dimensions
        ml_Rl, nl_Rl = size(R[k][2i-3]); #n_Rl = n_Rr
        ml_Rr, nl_Rr = size(R[k][2i-2]);
        m0_Rl, n0_Rl = size(R[k][2i-1]);
        m0_Rr, n0_Rr = size(R[k][2i]);
        mr_Rl, nr_Rl = size(R[k][2i+1]);
        mr_Rr, nr_Rr = size(R[k][2i+2]);

        ml_R̃l, nl_R̃l = size(R̃[k][2i-3]);
        ml_R̃r, nl_R̃r = size(R̃[k][2i-2]);
        m0_R̃l, n0_R̃l = size(R̃[k][2i-1]);
        m0_R̃r, n0_R̃r = size(R̃[k][2i]);
        mr_R̃l, nr_R̃l = size(R̃[k][2i+1]);
        mr_R̃r, nr_R̃r = size(R̃[k][2i+2]);

        ml_Wl, nl_Wl = size(W[k][2i-3]); #m_Wl = m_Wr
        ml_Wr, nl_Wr = size(W[k][2i-2]);
        m0_Wl, n0_Wl = size(W[k][2i-1]);
        m0_Wr, n0_Wr = size(W[k][2i]);
        mr_Wl, nr_Wl =size(W[k][2i+1]);
        mr_Wr, nr_Wr =size(W[k][2i+2]);

        ml_W̃l, nl_W̃l = size(W̃[k][2i-3]);
        ml_W̃r, nl_W̃r = size(W̃[k][2i-2]);
        m0_W̃l, n0_W̃l = size(W̃[k][2i-1]);
        m0_W̃r, n0_W̃r = size(W̃[k][2i]);
        mr_W̃l, nr_W̃l = size(W̃[k][2i+1]);
        mr_W̃r, nr_W̃r = size(W̃[k][2i+2]);

        mrr_Rl, nrr_Rl =size(R[k][2i+3]);
        mrr_Rr, nrr_Rr =size(R[k][2i+4]);
        mrr_R̃l, nrr_R̃l =size(R̃[k][2i+3]);
        mrr_R̃r, nrr_R̃r =size(R̃[k][2i+4]);
        mrr_Wl, nrr_Wl =size(W[k][2i+3]);
        mrr_Wr, nrr_Wr =size(W[k][2i+4]);
        mrr_W̃l, nrr_W̃l =size(W̃[k][2i+3]);
        mrr_W̃r, nrr_W̃r =size(W̃[k][2i+4]);

        mrrr_Wl, mrrr_Wr = size(W[k][2i+5])
        mrrr_R̃l, mrrr_R̃r = size(R̃[k][2i+5])

        #The size of the zero blocks in B correspond to the number of rows in each
        #coresponding block of R_C, and the number of columns in each block of W_C^H.

        #B_C[k][2i-1,2i+2]
        B_C[k][2i-1][2i+2] = [zeros(m0_Rl,mr_W̃r) B[k][2i-1][2i+1] B[k][2i-1][2i+2] zeros(m0_Rl,mrr_Wl);
                              zeros(ml_R̃r,mr_W̃r) zeros(ml_R̃r,mr_Wl) zeros(ml_R̃r,mr_Wr) zeros(ml_R̃r,mrr_Wl);
                              B̃[k][2i-1][2i+2] zeros(m0_R̃l,mr_Wl) zeros(m0_R̃l,mr_Wr) zeros(m0_R̃l,mrr_Wl);
                              B̃[k][2i][2i+2] zeros(m0_R̃r,mr_Wl) zeros(m0_R̃r,mr_Wr) zeros(m0_R̃r,mrr_Wl)] +
                                  R_C[k][2i-1]*F[k-1][i][i+1]*W_C[k][2i+2]'


        #B_C[k][2i+2,2i-1]
        B_C[k][2i+2][2i-1] = [zeros(mr_Rr,m0_W̃l) zeros(mr_Rr,ml_Wr) B[k][2i+2][2i-1] B[k][2i+2][2i];
                              B̃[k][2i+1][2i-1] zeros(mr_R̃l,ml_Wr) zeros(mr_R̃l,m0_Wl) zeros(mr_R̃l,m0_Wr);
                              B̃[k][2i+2][2i-1] zeros(mr_R̃r,ml_Wr) zeros(mr_R̃r,m0_Wl) zeros(mr_R̃r,m0_Wr);
                              zeros(mrr_R̃l,m0_W̃l) zeros(mrr_R̃l,ml_Wr) zeros(mrr_R̃l,m0_Wl) zeros(mrr_R̃l,m0_Wr)] +
                                  R_C[k][2i+2]*F[k-1][i+1][i]*W_C[k][2i-1]';

        #B_C[k][2i-1,2i+3]
        B_C[k][2i-1][2i+3] = [B[k][2i-1][2i+1]*G[k][2i+1]*B̃[k][2i+1][2i+3] B[k][2i-1][2i+2] zeros(m0_Rl,mrr_Wl) zeros(m0_Rl,mrr_Wr);
                              zeros(ml_R̃r,mrr_W̃l) zeros(ml_R̃r,mr_Wr) zeros(ml_R̃r,mrr_Wl) zeros(ml_R̃r,mrr_Wr);
                              zeros(m0_R̃l,mrr_W̃l) zeros(m0_R̃l,mr_Wr) zeros(m0_R̃l,mrr_Wl) zeros(m0_R̃l,mrr_Wr);
                              zeros(m0_R̃r,mrr_W̃l) zeros(m0_R̃r,mr_Wr) zeros(m0_R̃r,mrr_Wl) zeros(m0_R̃r,mrr_Wr)] +
                                  R_C[k][2i-1]*F[k-1][i][i+2]*W_C[k][2i+3]';

        #B_C[k][2i-1,2i+4]
        B_C[k][2i-1][2i+4] = [B[k][2i-1][2i+1]*G[k][2i+1]*B̃[k][2i+1][2i+4] + B[k][2i-1][2i+2]*G[k][2i+2]*B̃[k][2i+2][2i+4] zeros(m0_Rl,mrr_Wl) zeros(m0_Rl,mrr_Wr) zeros(m0_Rl,mrrr_Wl);
                              zeros(ml_R̃r,mrr_W̃r) zeros(ml_R̃r,mrr_Wl) zeros(ml_R̃r,mrr_Wr) zeros(ml_R̃r,mrrr_Wl);
                              zeros(m0_R̃l,mrr_W̃r) zeros(m0_R̃l,mrr_Wl) zeros(m0_R̃l,mrr_Wr) zeros(m0_R̃l,mrrr_Wl);
                              zeros(m0_R̃r,mrr_W̃r) zeros(m0_R̃l,mrr_Wl) zeros(m0_R̃l,mrr_Wr) zeros(m0_R̃l,mrrr_Wl)] +
                                  R_C[k][2i-1]*F[k-1][i][i+2]*W_C[k][2i+4]';

        #B_C[k][2i,2i+3]
        B_C[k][2i][2i+3] = [zeros(m0_Rr,mrr_W̃l) B[k][2i][2i+2] zeros(m0_Rr,mrr_Wl) zeros(m0_Rr,mrr_Wr);
                            zeros(m0_R̃l,mrr_W̃l) zeros(m0_R̃l,mr_Wr) zeros(m0_R̃l,mrr_Wl) zeros(m0_R̃l,mrr_Wr);
                            zeros(m0_R̃r,mrr_W̃l) zeros(m0_R̃r,mr_Wr) zeros(m0_R̃r,mrr_Wl) zeros(m0_R̃r,mrr_Wr);
                            B̃[k][2i+1][2i+3] zeros(mr_R̃l,mr_Wr) zeros(mr_R̃l,mrr_Wl) zeros(mr_R̃l,mrr_Wr)] +
                                R_C[k][2i]*F[k-1][i][i+2]*W_C[k][2i+3]';

        #B_C[k][2i,2i+4]
        B_C[k][2i][2i+4] = [B[k][2i][2i+2]*G[k][2i+2]*B̃[k][2i+2][2i+4] zeros(m0_Rr,mrr_Wl) zeros(m0_Rr,mrr_Wr) zeros(m0_Rr,mrrr_Wl);
                            zeros(m0_R̃l,mrr_W̃r) zeros(m0_R̃l,mrr_Wl) zeros(m0_R̃l,mrr_Wr) zeros(m0_R̃l,mrrr_Wl);
                            zeros(m0_R̃r,mrr_W̃r) zeros(m0_R̃r,mrr_Wl) zeros(m0_R̃r,mrr_Wr) zeros(m0_R̃r,mrrr_Wl);
                            B̃[k][2i+1][2i+4] zeros(mr_R̃l,mrr_Wl) zeros(mr_R̃l,mrr_Wr) zeros(mr_R̃l,mrrr_Wl)] +
                                R_C[k][2i]*F[k-1][i][i+2]*W_C[k][2i+4]';

#B_C[k][2i+3,2i-1]
B_C[k][2i+3][2i-1] = [B[k][2i+3][2i+1]*G[k][2i+1]*B̃[k][2i+1][2i-1] zeros(mrr_Rl,ml_Wr) zeros(mrr_Rl,m0_Wl) zeros(mrr_Rl,mr_Wl);
                      B̃[k][2i+2][2i-1] zeros(mr_R̃r,ml_Wr) zeros(mr_R̃r,m0_Wl) zeros(mr_R̃r,mr_Wl);
                      zeros(mrr_R̃l,m0_W̃l) zeros(mrr_R̃l,ml_Wr) zeros(mrr_R̃l,m0_Wl) zeros(mrr_R̃l,mr_Wl);
                      zeros(mrr_R̃r,m0_W̃l) zeros(mrr_R̃r,ml_Wr) zeros(mrr_R̃r,m0_Wl) zeros(mrr_R̃r,mr_Wl)] +
                          R_C[k][2i+3]*F[k-1][i+2][i]*W_C[k][2i-1]';

#B_C[k][2i+3,2i]
B_C[k][2i+3][2i] = [zeros(mrr_Rl,m0_W̃r) zeros(mrr_Rl,m0_Wl) zeros(mrr_Rl,m0_Wr) B[k][2i+3][2i+1];
                    B̃[k][2i+2][2i] zeros(mr_R̃r,m0_Wl) zeros(mr_R̃r,m0_Wr) zeros(mr_R̃r,mr_Wl);
                    zeros(mrr_R̃l,m0_W̃r) zeros(mrr_R̃l,m0_Wl) zeros(mrr_R̃l,m0_Wr) zeros(mrr_R̃l,mr_Wl);
                    zeros(mrr_R̃r,m0_W̃r) zeros(mrr_R̃r,m0_Wl) zeros(mrr_R̃r,m0_Wr) zeros(mrr_R̃r,mr_Wl)] +
                        R_C[k][2i+3]*F[k-1][i+2][i]*W_C[k][2i]';

#B_C[k][2i+4,2i-1]
B_C[k][2i+4][2i-1] = [B[k][2i+4][2i+1]*G[k][2i+1]*B̃[k][2i+1][2i-1] + B[k][2i+4][2i+2]*G[k][2i+2]*B̃[k][2i+2][2i-1] zeros(mrr_Rr,ml_Wr) zeros(mrr_Rr,m0_Wl) zeros(mrr_Rr,m0_Wr);
                      zeros(mrr_R̃l,m0_W̃l) zeros(mrr_R̃l,ml_Wr) zeros(mrr_R̃l,m0_Wl) zeros(mrr_R̃l,m0_Wr);
                      zeros(mrr_R̃r,m0_W̃l) zeros(mrr_R̃r,ml_Wr) zeros(mrr_R̃r,m0_Wl) zeros(mrr_R̃r,m0_Wr);
                      zeros(mrrr_R̃l,m0_W̃l) zeros(mrrr_R̃l,ml_Wr) zeros(mrrr_R̃l,m0_Wl) zeros(mrrr_R̃l,m0_Wr)] +
                          R_C[k][2i+4]*F[k-1][i+2][i]*W_C[k][2i-1]';

#B_C[k][2i+4,2i]
B_C[k][2i+4][2i] = [B[k][2i+4][2i+2]*G[k][2i+2]*B̃[k][2i+2][2i] zeros(mrr_Rr,m0_Wl) zeros(mrr_Rr,m0_Wr) B[k][2i+4][2i+1]
                    zeros(mrr_R̃l,m0_W̃r) zeros(mrr_R̃l,m0_Wl) zeros(mrr_R̃l,m0_Wr) zeros(mrr_R̃l,mr_Wl);
                    zeros(mrr_R̃r,m0_W̃r) zeros(mrr_R̃r,m0_Wl) zeros(mrr_R̃r,m0_Wr) zeros(mrr_R̃r,mr_Wl);
                    zeros(mrrr_R̃l,m0_W̃r) zeros(mrrr_R̃l,m0_Wl) zeros(mrrr_R̃l,m0_Wr) zeros(mrrr_R̃l,mr_Wl)] +
                        R_C[k][2i+4]*F[k-1][i+2][i]*W_C[k][2i]';

end
end

################################################################################
#Dense Blocks (D)

D_C = OffsetArray{OffsetArray{Array{Float64}}}(undef,1:2^depth+2);
for k = 1:2^depth+2
    D_C[k] = OffsetArray{Array{Float64}}(undef,k-2:k+2);
end

#Get sizes of F's to create empty U_C[depth^2+1] and V_C[depth^2+1] so the multiply will go through
k = depth

m_V1, n_V1 = size(F[k][2^depth-1][2^depth+1]);
m_U1, n_U1 = size(F[k][2^depth+1][2^depth-1]);
m_D̃1, n_D̃1 = size(D̃[2^depth][2^depth+1]);#size(D̃[k^2][k^2+1]);  these are wrong somehow -1?
m_D1, n_D1 = size(D[2^depth+1][2^depth]); #size(D[k^2+1][k^2]);
V_C[2^depth+1] = zeros(n_D̃1,n_V1); #U_C[k^2+1]
U_C[2^depth+1] = zeros(m_D1,m_U1); #U_C[k^2+1]

for i =1:2^depth-1 #k^2-1
    D_C[i][i] = D[i][i-1]*D̃[i-1][i] + D[i][i]*D̃[i][i] + D[i][i+1]D̃[i+1][i] + U_C[i]*F[k][i][i]*V_C[i]';
    D_C[i][i+1] = D[i][i]*D̃[i][i+1] + D[i][i+1]*D̃[i+1][i+1] + U_C[i]*F[k][i][i+1]*V_C[i+1]';
    D_C[i+1][i] = D[i+1][i]*D̃[i][i] + D[i+1][i+1]*D̃[i+1][i] + U_C[i+1]*F[k][i+1][i]*V_C[i]';
    D_C[i][i+2] = D[i][i+1]*D̃[i+1][i+2] + U_C[i]*F[k][i][i+2]*V_C[i+2]';
    D_C[i+2][i] = D[i+2][i+1]*D̃[i+1][i] + U_C[i+2]*F[k][i+2][i]*V_C[i]';
end

i = 2^depth; #k^2;
D_C[i][i] = D[i][i-1]*D̃[i-1][i] + D[i][i]*D̃[i][i] + D[i][i+1]D̃[i+1][i] + U_C[i]*F[k][i][i]*V_C[i]';

D_C,U_C,R_C,B_C,W_C,V_C

end


function Copy_FMM(fmm :: FMM)
    #Copy all matrices from the fmm structure into arrays
    #Input:             fmm: FMM
    #Output:            D:: OffsetArray{OffsetArray{Array{Float64}}}
    #                   U,V:: OffsetArray{Array{Float64}}
    #                   R,W:: Array{OffsetArray{Array{Float64}}}
    #                   B:: Array{OffsetArray{OffsetArray{Array{Float64}}}}

    ############################################################################
    #Copy U's V's
    U = OffsetArray{Array{Float64}}(undef,0:2^(fmm.depth)+2);
    V = OffsetArray{Array{Float64}}(undef,0:2^(fmm.depth)+2);

    i = 1;
    U,V = Copy_UV(fmm,U,V,i);

    #assign empty matrices to out of bounds U's, V's
    U[0] = zeros(0,0)
    U[2^(fmm.depth)+1] = zeros(0,0)
    V[0] = zeros(0,0)
    V[2^(fmm.depth)+1] = zeros(0,0)
    ############################################################################
    #Copy D's
    D = OffsetArray{OffsetArray{Array{Float64}}}(undef,0:2^fmm.depth+1);
    for k = 0:2^fmm.depth+1
        D[k] = OffsetArray{Array{Float64}}(undef,k-1:k+1); #(0:k+1)
    end
    i = 1;
    D = Copy_D(fmm,D,i,fmm.depth)

    #assign empty matrices to out of bounds D's
    D[1][0] = zeros(size(U[1],1),0);
    D[2^fmm.depth][2^fmm.depth+1] = zeros(size(U[2^fmm.depth],1),0);
    D[0][1] = zeros(0,size(V[1]',2));
    D[2^fmm.depth+1][2^fmm.depth] = zeros(0,size(V[2^fmm.depth]',2));

    ################################################################################
    #Copy R's, W's
    R = Array{OffsetArray{Array{Float64}}}(undef,fmm.depth);
    W = Array{OffsetArray{Array{Float64}}}(undef,fmm.depth);
    for k = 1:fmm.depth
        R[k] = OffsetArray{Array{Float64}}(undef,-1:2^(k)+4);
        W[k] = OffsetArray{Array{Float64}}(undef,-1:2^(k)+4);
    end
    i = 1;
    l = 1;
    R,W = Copy_RW(fmm,R,W,l,i)

    #assign empty R, W matrices to the right and left of the edge of each tree
    #(matrices outside the indices 1:2^k are 0x0 empty matrices)
    for k = 1:fmm.depth
        R[k][-1] = zeros(0,0);
        R[k][0] = zeros(0,0);
        R[k][2^k+1] = zeros(0,0);
        R[k][2^k+2] = zeros(0,0);
        R[k][2^k+3] = zeros(0,0);
        R[k][2^k+4] = zeros(0,0);
        W[k][-1] = zeros(0,0);
        W[k][0] = zeros(0,0);
        W[k][2^k+1] = zeros(0,0);
        W[k][2^k+2] = zeros(0,0);
        W[k][2^k+3] = zeros(0,0);
        W[k][2^k+4] = zeros(0,0);
    end
    ################################################################################
    #Copy B's

    B = Array{OffsetArray{OffsetArray{Array{Float64}}}}(undef,fmm.depth);
    for k = 1:fmm.depth
        B[k] = OffsetArray{OffsetArray{Array{Float64}}}(undef,-1:2^k+4);
        for i = -1:2^k+4
            if mod(i,2) == 1
                B[k][i] = OffsetArray{Array{Float64}}(undef,i-2:i+3); #(0:k+1)
            else
                B[k][i] = OffsetArray{Array{Float64}}(undef,i-3:i+2);
            end
        end
    end

    i = 1;
    k = 0;
    B = Copy_B(fmm,B,k,i)

    ################################################################################
    #assign empty B matrices to the right and left of the edge of each tree
    #(matrices outside the indices 1:2^k are ?x0 empty matrices)
    for k = 1:fmm.depth #1:fmm.depth-1
        k_1r, _ = size(R[k][1]); #size(R[k+1][1]);
        l_1r, _ = size(R[k][2]);
        k_1w, _ = size(W[k][1]); #~
        l_1w, _ = size(W[k][2]); #~

        k_3r, _ = size(R[k][2^(k)-1])
        l_3r, _ = size(R[k][2^(k)])
        k_3w, _ = size(W[k][2^(k)-1]) #~
        l_3w, _ = size(W[k][2^(k)]) #~


        B[k][-1][1] = zeros(0,k_1w);       B[k][1][-1] = zeros(k_1r,0);
        B[k][-1][2] = zeros(0,l_1w);       B[k][2][-1] = zeros(k_1r,0);
        B[k][0][2] = zeros(0,l_1w);        B[k][2][0] = zeros(l_1r,0);

        #print("k = ", k, "\n")
        B[k][2^(k)-1][2^(k)+1] = zeros(k_3r,0);  B[k][2^(k)+1][2^(k)-1] = zeros(0,k_3w);
        B[k][2^(k)-1][2^(k)+2] = zeros(k_3r,0);  B[k][2^(k)+2][2^(k)-1] = zeros(0,k_3w);
        B[k][2^(k)][2^(k)+2] = zeros(l_3r,0);    B[k][2^(k)+2][2^(k)] = zeros(0,l_3w);

        B[k][2^(k)+1][2^(k)+3] = zeros(0,0);     B[k][2^(k)+3][2^(k)+1] = zeros(0,0);
        B[k][2^(k)+1][2^(k)+4] = zeros(0,0);     B[k][2^(k)+4][2^(k)+1] = zeros(0,0);
        B[k][2^(k)+2][2^(k)+4] = zeros(0,0);     B[k][2^(k)+4][2^(k)+2] = zeros(0,0);
    end

    D,U,R,B,W,V

end

function FivePointFMMtoMatrix(D_C::OffsetArray{OffsetArray{Array{Float64}}},U_C::Array{Array{Float64}},R_C::Array{Array{Array{Float64}}},B_C::Array{OffsetArray{OffsetArray{Array{Float64}}}},W_C::Array{Array{Array{Float64}}},V_C::Array{Array{Float64}},fmm_idx::Array{Array{Int64}},depth::Int64)
    # Reconstruct the matrix C to test
    #Input:             D_C:: OffsetArray{OffsetArray{Array{Float64}}}
    #                   U_C,V_C:: Array{Array{Float64}}
    #                   R_C,W_C:: Array{Array{Array{Float64}}}
    #                   B_C:: Array{Array{Array{Float64}}}
    #                   fmm_idx:: Array{Array{Float64}} (stores index of each diagonal block at every level of the fmm structure)
    #                   depth:: Int64 (depth of the FMM structure)
    #Output:            C̃:: Array{Float64}


    BigU_C = Array{Array{Array{Float64}}}(undef,depth);
    BigV_C = Array{Array{Array{Float64}}}(undef,depth);


    for k = 1:depth
        BigU_C[k] = Array{Array{Float64}}(undef,2^k);
        BigV_C[k] = Array{Array{Float64}}(undef,2^k);
    end

    #Compute Big U's and V's
    for i = 1:2^depth
        BigU_C[depth][i] = U_C[i];
        BigV_C[depth][i] = V_C[i];
    end


    for k = depth:-1:2
        for i = 1:2^(k-1)
            BigU_C[k-1][i] = [BigU_C[k][2i-1]*R_C[k][2i-1];
                              BigU_C[k][2i]*R_C[k][2i]];
            BigV_C[k-1][i] = [BigV_C[k][2i-1]*W_C[k][2i-1];
                              BigV_C[k][2i]*W_C[k][2i]];
        end
    end

    #Multiply all the Blocks and insert them into their proper positions in C̃
    C̃ = zeros(N,N);

    # Add the last out of bounds index so the UBV' for loop below is easier to read
    for k = 1:depth+1
        fmm_idx[k] = [fmm_idx[k]; N+1];
    end

    #Insert the D blocks
    idx = fmm_idx[depth+1];

    for i = 1:2^depth

        if i == 1
            #first row has only 3 D blocks
            C̃[idx[1]:idx[2]-1,idx[1]:idx[2]-1] = D_C[1][1];
            C̃[idx[1]:idx[2]-1,idx[2]:idx[3]-1] = D_C[1][2];
            C̃[idx[1]:idx[2]-1,idx[3]:idx[4]-1] = D_C[1][3];

        elseif i == 2
            #second row has only 4 D blocks
            C̃[idx[2]:idx[3]-1,idx[1]:idx[2]-1] = D_C[2][1];
            C̃[idx[2]:idx[3]-1,idx[2]:idx[3]-1] = D_C[2][2];
            C̃[idx[2]:idx[3]-1,idx[3]:idx[4]-1] = D_C[2][3];
            C̃[idx[2]:idx[3]-1,idx[4]:idx[5]-1] = D_C[2][4];

        elseif i == 2^depth-1
            #second to last row has only 4 D blocks
            C̃[idx[i]:idx[i+1]-1,idx[i-2]:idx[i-1]-1] = D_C[i][i-2];
            C̃[idx[i]:idx[i+1]-1,idx[i-1]:idx[i]-1] = D_C[i][i-1];
            C̃[idx[i]:idx[i+1]-1,idx[i]:idx[i+1]-1] = D_C[i][i];
            C̃[idx[i]:idx[i+1]-1,idx[i+1]:idx[i+2]-1] = D_C[i][i+1];

        elseif i == 2^depth
            #last row has only 3 D blocks
            C̃[idx[i]:idx[i+1]-1,idx[i-2]:idx[i-1]-1] = D_C[i][i-2];
            C̃[idx[i]:idx[i+1]-1,idx[i-1]:idx[i]-1] = D_C[i][i-1];
            C̃[idx[i]:idx[i+1]-1,idx[i]:idx[i+1]-1] = D_C[i][i];

        else
            #all other rows have 5 D blocks
            C̃[idx[i]:idx[i+1]-1,idx[i-2]:idx[i-1]-1] = D_C[i][i-2];
            C̃[idx[i]:idx[i+1]-1,idx[i-1]:idx[i]-1] = D_C[i][i-1];
            C̃[idx[i]:idx[i+1]-1,idx[i]:idx[i+1]-1] = D_C[i][i];
            C̃[idx[i]:idx[i+1]-1,idx[i+1]:idx[i+2]-1] = D_C[i][i+1];
            C̃[idx[i]:idx[i+1]-1,idx[i+2]:idx[i+3]-1] = D_C[i][i+2];

        end

    end

    #Multiply and insert all UBV' blocks

    #for each level
    for k = 2:size(B_C,1)
        #step through each L shaped block
        for i = 1:2^(k-1)-1 #-1 because nothing to do at the last block

            if i == 2^(k-1)-1
                #only 2 UVB blocks at the last L shaped blocks for every level
                C̃[fmm_idx[k+1][2i-1]:fmm_idx[k+1][2i]-1,fmm_idx[k+1][2i+2]:fmm_idx[k+1][2i+3]-1] = BigU_C[k][2i-1]*B_C[k][2i-1][2i+2]*BigV_C[k][2i+2]';

                C̃[fmm_idx[k+1][2i+2]:fmm_idx[k+1][2i+3]-1,fmm_idx[k+1][2i-1]:fmm_idx[k+1][2i]-1] = BigU_C[k][2i+2]*B_C[k][2i+2][2i-1]*BigV_C[k][2i-1]';

            else
                #10 UVB blocks at every other level
                C̃[fmm_idx[k+1][2i-1]:fmm_idx[k+1][2i]-1,fmm_idx[k+1][2i+2]:fmm_idx[k+1][2i+3]-1] = BigU_C[k][2i-1]*B_C[k][2i-1][2i+2]*BigV_C[k][2i+2]';
                C̃[fmm_idx[k+1][2i-1]:fmm_idx[k+1][2i]-1,fmm_idx[k+1][2i+3]:fmm_idx[k+1][2i+4]-1] = BigU_C[k][2i-1]*B_C[k][2i-1][2i+3]*BigV_C[k][2i+3]';
                C̃[fmm_idx[k+1][2i-1]:fmm_idx[k+1][2i]-1,fmm_idx[k+1][2i+4]:fmm_idx[k+1][2i+5]-1] = BigU_C[k][2i-1]*B_C[k][2i-1][2i+4]*BigV_C[k][2i+4]';
                C̃[fmm_idx[k+1][2i]:fmm_idx[k+1][2i+1]-1,fmm_idx[k+1][2i+3]:fmm_idx[k+1][2i+4]-1] = BigU_C[k][2i]*B_C[k][2i][2i+3]*BigV_C[k][2i+3]'; #here
                C̃[fmm_idx[k+1][2i]:fmm_idx[k+1][2i+1]-1,fmm_idx[k+1][2i+4]:fmm_idx[k+1][2i+5]-1] = BigU_C[k][2i]*B_C[k][2i][2i+4]*BigV_C[k][2i+4]';

                C̃[fmm_idx[k+1][2i+2]:fmm_idx[k+1][2i+3]-1,fmm_idx[k+1][2i-1]:fmm_idx[k+1][2i]-1] = BigU_C[k][2i+2]*B_C[k][2i+2][2i-1]*BigV_C[k][2i-1]';
                C̃[fmm_idx[k+1][2i+3]:fmm_idx[k+1][2i+4]-1,fmm_idx[k+1][2i-1]:fmm_idx[k+1][2i]-1] = BigU_C[k][2i+3]*B_C[k][2i+3][2i-1]*BigV_C[k][2i-1]';
                C̃[fmm_idx[k+1][2i+4]:fmm_idx[k+1][2i+5]-1,fmm_idx[k+1][2i-1]:fmm_idx[k+1][2i]-1] = BigU_C[k][2i+4]*B_C[k][2i+4][2i-1]*BigV_C[k][2i-1]'; #here
                C̃[fmm_idx[k+1][2i+3]:fmm_idx[k+1][2i+4]-1,fmm_idx[k+1][2i]:fmm_idx[k+1][2i+1]-1] = BigU_C[k][2i+3]*B_C[k][2i+3][2i]*BigV_C[k][2i]';
                C̃[fmm_idx[k+1][2i+4]:fmm_idx[k+1][2i+5]-1,fmm_idx[k+1][2i]:fmm_idx[k+1][2i+1]-1] = BigU_C[k][2i+4]*B_C[k][2i+4][2i]*BigV_C[k][2i]';
            end

        end

    end

C̃
end
################################################################################

function Copy_UV(fmm :: FMM, U :: OffsetArray{Array{Float64}}, V :: OffsetArray{Array{Float64}}, i :: Int64)
    if isa(fmm, Leaf)
        U[i] = fmm.U;
        V[i] = fmm.V;

        U,V

    elseif isa(fmm,Node)
        U1,V1 = Copy_UV(fmm.fmmUL,U,V,2*i-1)
        U2,V2 = Copy_UV(fmm.fmmLR,U1,V1,2*i)

        U2,V2

    end
end

function Copy_RW(fmm :: FMM, R :: Array{OffsetArray{Array{Float64}}}, W :: Array{OffsetArray{Array{Float64}}}, l :: Int64, i :: Int64)
    if isa(fmm.fmmUL, Leaf) #assuming a complete tree
        R[l][2i-1] = fmm.Rl;
        R[l][2i] = fmm.Rr;

        W[l][2i-1] = fmm.Wl;
        W[l][2i] = fmm.Wr;

        R,W
    elseif isa(fmm.fmmUL, Node) #assuming a complete tree
        R[l][2i-1] = fmm.Rl;
        R[l][2i] = fmm.Rr;

        W[l][2i-1] = fmm.Wl;
        W[l][2i] = fmm.Wr;
        #println("2i-1 = ", 2i-1)

        R1,W1 = Copy_RW(fmm.fmmUL,R,W,l+1,2i-1);
        R2,W2 = Copy_RW(fmm.fmmLR,R,W,l+1,2i);

        R2,W2
    end
end


function Copy_D(fmm :: FMM, D :: OffsetArray{OffsetArray{Array{Float64}}}, i :: Int64, depth)
    if isa(fmm, Leaf)

        if i == 1
            D[i][i] = fmm.D_0;
            D[i][i+1] = fmm.D_r;
        elseif i == 2 .^depth
            D[i][i-1] = fmm.D_l;
            D[i][i] = fmm.D_0;
        else
            D[i][i-1] = fmm.D_l;
            D[i][i] = fmm.D_0;
            D[i][i+1] = fmm.D_r;
        end

        D

    elseif isa(fmm,Node)
        D1 = Copy_D(fmm.fmmUL,D,2*i-1,depth)
        D2 = Copy_D(fmm.fmmLR,D1,2*i,depth)

        D2
    end
end

function Copy_B(fmm :: FMM, B :: Array{OffsetArray{OffsetArray{Array{Float64}}}}, k :: Int64, i :: Int64) #Copy_B(fmm :: FMM, B :: Array{OffsetArray{Array{Array{Float64}}}}, k :: Int64, i :: Int64)

    if isa(fmm.fmmUL, Leaf) #assuming a complete tree
        # No B's store at the leaf.  Return
        B

    elseif isa(fmm.fmmUR, Node) #assuming a complete tree

        B1 = Copy_B_OffDiag(fmm.fmmUR,fmm.fmmLL,B,k+1,i)
        B2 = Copy_B(fmm.fmmUL,B1,k+1,2i-1);
        B3 = Copy_B(fmm.fmmLR,B2,k+1,2(i+1)-1);

        B3
    end
end

function Copy_B_OffDiag(fmmUR :: FMM, fmmLL :: FMM, B :: Array{OffsetArray{OffsetArray{Array{Float64}}}}, k :: Int64, i :: Int64)#Copy_B_OffDiag(fmmUR :: FMM, fmmLL :: FMM, B :: Array{OffsetArray{Array{Array{Float64}}}}, k :: Int64, i :: Int64)
    if isa(fmmUR, LeafOffDiag) #assuming a complete tree

        #No B's stored at the leaf.  Return.
        B

    elseif isa(fmmUR, Node) #assuming a complete tree
        #this will place all B's between cousins in the array

        i_new = 2i-1;

        B[k+1][i_new][i_new+2] = fmmUR.B13;
        B[k+1][i_new][i_new+3] = fmmUR.B14;
        B[k+1][i_new+1][i_new+3] = fmmUR.B24;

        B[k+1][i_new+2][i_new] = fmmLL.B31;
        B[k+1][i_new+3][i_new] = fmmLL.B41;
        B[k+1][i_new+3][i_new+1] = fmmLL.B42;

        B1 = Copy_B_OffDiag(fmmUR.fmmLL,fmmLL.fmmUR,B,k+1,i_new+1)

        B1
    end

end

function get_fmm_indices(fmm :: FMM,fmm_idx :: Array{Array{Int64}},fmm_m :: Array{Array{Int64}},k :: Int64,i :: Int64, depth :: Int64)
#stores indices given in the fmm structure into an array of arrays
#INPUT:                fmm (FMM) fmm structure
#                      fmm_idx (Array{Array{Array{Int64}}}) Empty array of arrays
#                      k (Int64) level we are at in the tree/fmm structure
#                      i (Int64) current node at the kth level
#                      depth (Int64) total depth of the fmm tree
#OUTPUT                fmm_idx (Array{Array{Int64}}) Array of arrays containing fmm indices
#                      fmm_idx (Array{Array{Int64}}) Array of arrays containing fmm partition dimensions

    if isa(fmm, Leaf)

        fmm_idx[depth-k+1][i] = fmm.col_idx;
        fmm_m[depth-k+1][i] = fmm.m;

        fmm_idx,fmm_m
    elseif isa(fmm, Node)

        fmm_idx[depth-k+2][2i-1] = fmm.fmmUL.col_idx;
        fmm_m[depth-k+2][2i-1] = fmm.fmmUL.m;

        fmm_idx[depth-k+2][2i] = fmm.fmmLR.col_idx;
        fmm_m[depth-k+2][2i] = fmm.fmmLR.m;

        fmm_idx,fmm_m = get_fmm_indices(fmm.fmmUL,fmm_idx,fmm_m,k-1,2i-1,depth)
        fmm_idx,fmm_m = get_fmm_indices(fmm.fmmLR,fmm_idx,fmm_m,k-1,2i,depth)

        fmm_idx,fmm_m
    end

end
