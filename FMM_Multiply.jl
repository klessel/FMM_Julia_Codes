function FMM_Multiply(fmm::FMM,x::Array{Float64,2}) #(fmm::FMM,x::Array{Float64}(2))
#function [b] = HSS_Multiply(hss,x)
#Takes in an a matrix in fmm form, A_fmm, and multiplies this with an input
#array, x.  Gives an array b as output.  A_fmm*x = b.
#
#INPUT:             fmm     (FMM) a matrix in fmm form as given by the
#                           function FMM_Construction(). Contains partition
#                           dimensions and matrices (U, V, R,W and B), which are
#                           heirachically stored, and compose the fmm
#                           representation of the input matrix.
#                   x       (Array{Float64,2}) input array to multiply with A.
#
#OUTPUT:            b       (Array{Float64,2}) output array from the resulting multiplication
#
# Author: Kristen Lessel, klessel@engineering.ucsb.edu
#
# Copyright (C) 2016 Kristen Lessel, (klessel@engineering.ucsb.edu)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

    # Notes regarding notation:
    # fmm0 contains 4 substructures, fmmUL, fmmUR, fmmLL and fmmLR, shown here:
    #     fmm0
    #   ---------
    #  \ UL \ UR \
    #  \----\----\
    #  \ LL \ LR \
    #   ---------
    #
    # Visual interpretation of a left call:
    #    fmmL            fmm0           fmmR
    #   -------        -------        -------
    #  \   \ x \      \ x \ x \      \   \   \
    #  \---\---\      \---\---\      \---\---\
    #  \   \   \      \   \   \      \   \   \
    #   -------        -------        -------
    # Next level:
    #            fmmL            fmm0   fmmR
    #   ---------------        ---------------
    #  \   \   \   \ x \      \ x \ x \   \   \
    #  \-------\-------\      \-------\-------\
    #  \   \   \   \   \      \   \   \   \   \
    #  \---------------\      \---------------\
    #  \   \   \   \   \      \   \   \   \   \
    #  \---------------\      \-------\-------\
    #  \   \   \   \   \      \   \   \   \   \
    #   ---------------        ---------------
    # So when going down one level in the fmm tree we set our new fmm0 = fmm0.fmmUL,
    # set the new fmmL =  fmmL.fmmUR, and set the new fmmR = fmm0.fmmUR
    # (right calls work similarly)
    # Visual interpretation of a right call:
    #    fmmL            fmm0           fmmR
    #   -------        -------        -------
    #  \   \   \      \   \   \      \   \   \
    #  \---\---\      \---\---\      \---\---\
    #  \   \   \      \ x \ x \      \ x \   \
    #   -------        -------        -------
    # Next level:
    #     fmmL   fmmO            fmmR
    #   ---------------        ---------------
    #  \   \   \   \   \      \   \   \   \   \
    #  \-------\-------\      \-------\-------\
    #  \   \   \   \   \      \   \   \   \   \
    #  \---------------\      \---------------\
    #  \   \   \   \   \      \   \   \   \   \
    #  \-------\-------\      \-------\-------\
    #  \   \   \ x \ x \      \ x \   \   \   \
    #   ---------------        ---------------


    f = Array{Float64}(undef,0,1); # Array(Float64, size(fmm.Rl,2),1);

    gtree = FMM_Upsweep(fmm,x);

    depth = 0;
    empty_leafoffdiag = LeafOffDiag(-1,-1,-1,-1,depth,Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                             Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),Array{Float64}(undef,0,0),
                             Array{Float64}(undef,0,0));

    g = Array{Float64}(undef,0,1);
    empty_gleaf = Leaf_Gtree(0,1,depth,g);

    b = FMM_Downsweep(empty_leafoffdiag,fmm,empty_leafoffdiag,x,f,empty_gleaf,gtree,empty_gleaf);

    b

end

function FMM_Upsweep(fmm,x)
#upsweep to compute g's

if isa(fmm,Leaf)

    g = fmm.V'*x;
    gtree = Leaf_Gtree(fmm.m,fmm.col_idx,fmm.depth,g);

    gtree
else
    xl = x[1:fmm.fmmUL.m,:];
    gtreeL = FMM_Upsweep(fmm.fmmUL,xl);

    xr = x[fmm.fmmUL.m+1:end,:];
    gtreeR = FMM_Upsweep(fmm.fmmLR,xr);

    #new g
    g = fmm.Wl'*gtreeL.g+fmm.Wr'*gtreeR.g;

    gtree = Node_Gtree(gtreeL,gtreeR,fmm.m,fmm.col_idx,fmm.depth,g)

    gtree
end

end

function FMM_Downsweep(fmmL,fmm0,fmmR,x,f,gtreeL,gtree0,gtreeR)
#downsweep to compute f's

if isa(gtreeL, Leaf_Gtree) && isa(gtree0, Leaf_Gtree) && isa(gtreeR, Leaf_Gtree)
# Case 1) if all are leaves

    x_l = x[gtreeL.col_idx:(gtreeL.col_idx+gtreeL.m-1),:];
    x_0 = x[gtree0.col_idx:(gtree0.col_idx+gtree0.m-1),:];
    x_r = x[gtreeR.col_idx:(gtreeR.col_idx+gtreeR.m-1),:];

    b = fmm0.U*f + fmm0.D_l*x_l + fmm0.D_0*x_0 + fmm0.D_r*x_r;
    b

#elseif isa(fmm0, Node)
elseif isa(gtreeL, Node_Gtree) && isa(gtree0,Node_Gtree) && isa(gtreeR,Node_Gtree)
    # Case 2) if all are nodes
    #      fmmL    fmm0    fmmR
    #     -------------------------------
    #    \ . \ . \   \   \       \       \
    #    \-------\-------\       \       \
    #    \ . \ . \ . \   \       \       \
    #     ---------------\---------------
    #    \B31\ . \ . \ . \B13\B14\       \
    # 2) \-------\-------\-------\       \
    #    \B41\B42\ . \ . \ . \B24\       \
    #     -------------------------------
    #    \       \   \ . \ . \ . \   \   \
    #    \       \---\---\---\---\---\---\
    #    \       \   \   \ . \ . \ . \   \
    #     ---------------\---------------
    #    \       \       \   \ . \ . \ . \
    #    \       \       \-------\-------\
    #    \       \       \   \   \ . \ . \
    #     -------------------------------

    fl =  fmm0.Rl*f + fmmL.B31*gtreeL.gtreeL.g + fmmR.B13*gtreeR.gtreeL.g + fmmR.B14*gtreeR.gtreeR.g;
    bl = FMM_Downsweep(fmmL.fmmUR,fmm0.fmmUL,fmm0.fmmUR,x,fl,gtreeL.gtreeR,gtree0.gtreeL,gtree0.gtreeR);

    fr = fmm0.Rr*f + fmmL.B41*gtreeL.gtreeL.g + fmmL.B42*gtreeL.gtreeR.g + fmmR.B24*gtreeR.gtreeR.g;
    br = FMM_Downsweep(fmm0.fmmLL,fmm0.fmmLR,fmmR.fmmLL,x,fr,gtree0.gtreeL,gtree0.gtreeR,gtreeR.gtreeL)
    b = [bl; br];

    b

elseif isa(gtreeL, Node_Gtree) && isa(gtree0, Node_Gtree) && isa(gtreeR, Leaf_Gtree)
    #Case 3) if Left and Center Trees are nodes, Right is Leaf
    #      fmmL    fmm0    fmmR
    #     -------------------------------
    #    \ . \ . \   \   \       \       \
    #    \-------\-------\       \       \
    #    \ . \ . \ . \   \       \       \
    #     ---------------\---------------
    #    \   \ . \ . \ . \  B13  \       \
    # 3) \-------\-------\-------\       \
    #    \   \   \ . \ . \   .   \       \
    #     -------------------------------
    #    \       \   \   \       \   \   \
    #    \       \   \ . \   .   \ . \   \
    #    \       \   \   \       \   \   \
    #     ---------------\---------------
    #    \       \       \   \ . \ . \ . \
    #    \       \       \-------\-------\
    #    \       \       \   \   \ . \ . \
    #     -------------------------------

    #fix addition of 0x0 (empty) arrays to px1 arrays (for the second call of the routine each B*g needs to be of size px1 for the addition to go through)
    if isempty(fmmR.B13)
        #If we are at the right most branch of the tree
        fl = fmm0.Rl*f + fmmL.B31*gtreeL.gtreeL.g;
    else
        fl = fmm0.Rl*f + fmmL.B31*gtreeL.gtreeL.g + fmmR.B13*gtreeR.g;
    end
    bl = FMM_Downsweep(fmmL.fmmUR,fmm0.fmmUL, fmm0.fmmUR,x,fl,gtreeL.gtreeR,gtree0.gtreeL,gtree0.gtreeR);

    fr = fmm0.Rr*f + fmmL.B41*gtreeL.gtreeL.g + fmmL.B42*gtreeL.gtreeR.g;
    br = FMM_Downsweep(fmm0.fmmLL,fmm0.fmmLR,fmmR,x,fr,gtree0.gtreeL,gtree0.gtreeR,gtreeR);

    b = [bl;br];

elseif isa(gtreeL, Node_Gtree) && isa(gtree0, Leaf_Gtree) && isa(gtreeR, Node_Gtree)
    #Case 4) if Left and Right Trees are nodes, Center is Leaf
    #              fmmL     fmm0   fmmR
    #     -------------------------------
    #    \ . \ . \   \   \       \       \
    #    \-------\-------\       \       \
    #    \ . \ . \ . \   \       \       \
    #     ---------------\---------------
    #    \   \ . \ . \ . \       \       \
    #    \-------\-------\-------\       \
    #    \   \   \ . \ . \   .   \       \
    #     -------------------------------
    #    \       \   \   \       \   \   \
    # 4) \       \B31\ . \   .   \ . \B13\
    #    \       \   \   \       \   \   \
    #     ---------------\---------------
    #    \       \       \   .   \ . \ . \
    #    \       \       \-------\-------\
    #    \       \       \       \ . \ . \
    #     -------------------------------

    # Here R = I
    f_new = f + fmmL.B31*gtreeL.gtreeL.g + fmmR.B13*gtreeR.gtreeR.g;
    b = FMM_Downsweep(fmmL.fmmUR,fmm0,fmmR.fmmLL,x,f_new,gtreeL.gtreeR,gtree0,gtreeR.gtreeL);

    b

elseif isa(gtreeL, Leaf_Gtree) && isa(gtree0,Node_Gtree) && isa(gtreeR,Node_Gtree)
    #Case 5) if left tree is a leaf, and center and right trees are nodes
    #              fmmL     fmm0   fmmR
    #     -------------------------------
    #    \       \       \       \       \
    #    \   .   \   .   \       \       \
    #    \       \       \       \       \
    #     ---------------\---------------
    #    \       \       \   \   \       \
    #    \   .   \   .   \ . \   \       \
    #    \       \       \   \   \       \
    #     -------------------------------
    #    \       \   .   \ . \ . \   \   \
    # 5) \       \-------\-------\-------\
    #    \       \  B31  \ . \ . \ . \   \
    #     ---------------\---------------
    #    \       \       \   \ . \ . \ . \
    #    \       \       \-------\-------\
    #    \       \       \   \   \ . \ . \
    #     -------------------------------

    fl = fmm0.Rl*f + fmmR.B13*gtreeR.gtreeL.g + fmmR.B14*gtreeR.gtreeR.g;
    bl = FMM_Downsweep(fmmL,fmm0.fmmUL,fmm0.fmmUR,x,fl,gtreeL,gtree0.gtreeL,gtree0.gtreeR);

    #fix addition of 0x0 arrays to px1 arrays (for the second call of the routine each B*g needs to be of size px1 for the addition to go through)
    if isempty(fmmL.B31)
        #If we are at the leftmost branch of the tree
        fr = fmm0.Rr*f + fmmR.B24*gtreeR.gtreeR.g;
    else
        fr = fmm0.Rr*f + fmmL.B31*gtreeL.g + fmmR.B24*gtreeR.gtreeR.g;
    end
    br = FMM_Downsweep(fmm0.fmmLL,fmm0.fmmLR,fmmR.fmmLL,x,fr,gtree0.gtreeL,gtree0.gtreeR,gtreeR.gtreeL);

    b = [bl;br];

    b

elseif isa(gtreeL, Leaf_Gtree) && isa(gtree0,Leaf_Gtree) && isa(gtreeR,Node_Gtree)
    #Case 6) if left and center trees are leaves, and right tree is a node
    #       fmmL   fmm0    fmmR
    #     -------------------------------
    #    \       \       \       \       \
    #    \   .   \   .   \       \       \
    #    \       \       \       \       \
    #     ---------------\---------------
    #    \       \       \   \   \       \
    # 6) \   .   \   .   \ . \B13\       \
    #    \       \       \   \   \       \
    #     -------------------------------
    #    \       \   .   \ . \ . \   \   \
    #    \       \-------\-------\-------\
    #    \       \       \ . \ . \ . \   \
    #     ---------------\---------------
    #    \       \       \   \ . \ . \ . \
    #    \       \       \-------\-------\
    #    \       \       \   \   \ . \ . \
    #     -------------------------------

    # Here R = I
    f_new = f + fmmR.B13*gtreeR.gtreeR.g;
    b = FMM_Downsweep(fmmL,fmm0,fmmR.fmmLL,x,f_new,gtreeL,gtree0,gtreeR.gtreeL);

    b

elseif isa(gtreeL, Node_Gtree) && isa(gtree0,Leaf_Gtree) && isa(gtreeR,Leaf_Gtree)
    #Case 7) if left tree is a Node, and center and right trees are leaves
    #       fmmL   fmm0    fmmR
    #     -------------------------------
    #    \ . \ . \       \       \       \
    #    \-------\-------\       \       \
    #    \ . \ . \   .   \       \       \
    #     ---------------\---------------
    #    \   \   \       \       \       \
    # 7) \B31\ . \   .   \   .   \       \
    #    \   \   \       \       \       \
    #     -------------------------------
    #    \       \       \       \   \   \
    #    \       \   .   \   .   \ . \   \
    #    \       \       \       \   \   \
    #     ---------------\---------------
    #    \       \       \   .   \ . \ . \
    #    \       \       \-------\-------\
    #    \       \       \       \ . \ . \
    #     -------------------------------

    f_new = f + fmmL.B31*gtreeL.gtreeL.g;
    b = FMM_Downsweep(fmmL.fmmUR,fmm0,fmmR,x,f_new,gtreeL.gtreeR,gtree0,gtreeR);

    b

elseif isa(gtreeL, Leaf_Gtree) && isa(gtree0,Node_Gtree) && isa(gtreeR,Leaf_Gtree)
    #Case 8) if left and right trees are Nodes, and center tree is a leaf
    # At the root both the left and the right leaves will be emtpy
    #       fmmL   fmm0    fmmR
    #     -------------------------------
    #    \       \   \   \       \       \
    #    \   .   \ . \   \       \       \
    #    \       \   \   \       \       \
    #     ---------------\---------------
    #    \   .   \ . \ . \  B13  \       \
    # 8) \-------\-------\-------\       \
    #    \  B31  \ . \ . \   .   \       \
    #     -------------------------------
    #    \       \   \   \       \       \
    #    \       \   \ . \   .   \   .   \
    #    \       \   \   \       \       \
    #     ---------------\---------------
    #    \       \       \       \       \
    #    \       \       \   .   \   .   \
    #    \       \       \       \       \
    #     -------------------------------

    if isempty(fmmR.B13)
        #We are at the right most branch of the tree ( at the root of the tree
        #these dimensions match, but with an uneven tree we have to modify the equation
        #for empty LeafOffDiags)
        fl = fmm0.Rl*f;
    else
        fl = fmm0.Rl*f + fmmR.B13*gtreeR.g;
    end
    bl = FMM_Downsweep(fmmL,fmm0.fmmUL,fmm0.fmmUR,x,fl,gtreeL,gtree0.gtreeL,gtree0.gtreeR);

    if isempty(fmmL.B31)
       #if we are at the leftmost branch of the tree ( at the root of the tree
        #these dimensions match, but with an uneven tree we have to modify the equation
        #for empty LeafOffDiags)
        fr = fmm0.Rr*f;
    else
        fr = fmm0.Rr*f + fmmL.B31*gtreeL.g;
    end
    br = FMM_Downsweep(fmm0.fmmLL,fmm0.fmmLR,fmmR,x,fr,gtree0.gtreeL,gtree0.gtreeR,gtreeR);

    b = [bl;br];
    b
end

end
