# Display contents of the packages since
# we might use them latter
with(LinearAlgebra):
with(Groebner):

writeto("output.txt"):

#F := [x^3 + 12, 2*y^2 - x]:
#a, b := Basis(F, plex(x, y), output=extended):

#simplify(Multiply(convert(b, Matrix), Transpose(convert(F, Matrix)))):
#simplify(Transpose(convert(a, Matrix))):

#_m_order_1 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]:
#order_1 = 'matrix'(_m_order_1, [x, y, z]):


#f_example_monomial_order := proc (poly, order_1::{MonomialOrder, ShortMonomialOrder})
#return sort([op(poly)], (b, a) -> TestOrder(a, b, order_1))
#end proc:

#F := 1/2*x + 2*y*z:

#print(f_example_monomial_order(F, 'matrix'([[1,0,0],[0,1,0],[0,0,1]], [x, y, z])        )):
#print(f_example_monomial_order(F, plex(x, y, z)   )):
#print(f_example_monomial_order(F, grlex(x, y, z)  )):

# ---------------------------------------------------------------------------
truncatePolynomial := proc (poly, order_1, order_2)
local leading_coeff_1, leading_mon_1, leading_term, curr_index, polys:

    leading_coeff_1, leading_mon_1 := LeadingTerm(poly, order_1):
    leading_term                   := leading_coeff_1 * leading_mon_1:
    curr_index                     := 1:

    polys := sort([op(poly)], (b, a) -> TestOrder(a, b, order_2)):

    while evalb(polys[curr_index] <> leading_term) do
        curr_index := curr_index + 1:
    end do:

    if curr_index = 1 then
        return polys[1]:
    else
        return add(x, x in polys[1..curr_index]):
    end if:
end proc:
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Test truncatePolynomial
print( truncatePolynomial(x^2 - y*z + 1, grlex(x, y, z), plex(x, y, z)) ):
print( truncatePolynomial(y^2 - x, grlex(x, y, z), plex(x, y, z)) ):
# ---------------------------------------------------------------------------

basisConversion := proc (basis, order_1, order_2)
local F, F_t, num_iter, F_t_gb, multipliers, G, repeat, curr_index, num_elements:

    print("Step 1"):
    F   := Basis(basis, order_1):
    F_t := map(v -> truncatePolynomial(v, order_1, order_2), F):
    print("Current F", F):
    print("Current F_t", F_t):

    num_iter := 1:

    while true do
        print("Iteration ", num_iter):
        num_iter := num_iter + 1:

        print("Step 3 & 4"):
        F_t_gb, multipliers := Basis(F_t, order_2, output=extended):
        print("Matrix of multipliers M': ", multipliers):
        print("H (= M' F_t): ", F_t_gb):

        print("Step 5"):
        G := convert(simplify(Multiply(convert(multipliers, Matrix), Transpose(convert(F, Matrix)))), list):
        print("G (= M' * F): ", G):

        print("Step 6 and 7"):
        repeat := false:
        curr_index := 1:
        num_elements := numelems(F_t_gb):
        while curr_index < num_elements do
            if evalb(LeadingMonomial(F_t_gb[curr_index], order_2) <> LeadingMonomial(G[curr_index], order_2)) then
                print("There are the witness polynomials that prove G_s is not empty (g_j, h_j) ", F_t_gb[curr_index], G[curr_index]):
                repeat := true:
                break:
            end if:
            curr_index := curr_index + 1:
        end do:

        if repeat = false then
            print("Done"):
            return InterReduce(G, order_2):
        else
            F := G:
            F_t := map(v -> truncatePolynomial(v, order_1, order_2), F):
        end if:
    end do:

end proc:

# ---------------------------------------------------------------------------
# Test basisConversion
print( basisConversion([y^2-x,x^2-y*z-1,z^2-x], grlex(x, y, z), plex(x, y, z)) ):
# ---------------------------------------------------------------------------

