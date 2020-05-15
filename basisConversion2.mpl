with(LinearAlgebra, 
  Multiply, Transpose):
with(Groebner,      
  Basis, InterReduce, Reduce, LeadingTerm, LeadingMonomial, TestOrder): 

writeto("output.txt"):

truncatePolynomial := proc (poly, order_1, order_2)
local leading_coeff_1, leading_mon_1, leading_term, 
curr_index, polys:
 
if evalb(type(poly, `+`)) then
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
else
  return poly:
end if:

end proc:

basisConversion := proc (basis, order_1, order_2)
local F, F_t, aux_F_t, num_iter, 
F_t_gb, H_t, multipliers, 
G, repeat, curr_index, num_elements,
tru_g_i:

lprint("Step 1"):
F   := Basis(basis, order_1):
F_t := map(v -> truncatePolynomial(v, order_1, order_2), F):
lprint("Current F", F):
lprint("Current F_t", F_t):

num_iter := 1:

while true do
  lprint("Iteration ", num_iter):
  num_iter := num_iter + 1:

  lprint("Step 2"):
  F_t_gb, multipliers := Basis(F_t, order_2, output=extended):
  lprint("Matrix of multipliers M': ", multipliers):
  lprint("H (= M' F_t): ", F_t_gb):
  
  lprint("Step 3"):
  H_t := map(v -> truncatePolynomial(v, order_1, order_2), F_t_gb):
  lprint("H_t ", H_t):

  lprint("Step 4"):
  G := convert(simplify(
  Multiply(convert(multipliers, Matrix), Transpose(convert(F, Matrix)))), list):
  lprint("G (= M' * F): ", G):

  lprint("Step 5"):
  repeat       := false:
  curr_index   := 1:
  num_elements := numelems(H_t):
  while curr_index < num_elements do
    tru_g_i := truncatePolynomial(G[curr_index], order_1, order_2):
    if evalb(tru_g_i <> H_t[curr_index]) then
      repeat := true:
    end if:
    if evalb(G[curr_index] <> 0) then
      H_t := subsop(curr_index = tru_g_i, H_t):
    end if:
    curr_index := curr_index + 1;
  end do:

  #F := InterReduce(G, order_2):
  F := G:
  if repeat = false then
    lprint("Done"):
    #return F:
    return InterReduce(F, order_2):
  else
    #F_t := map(v -> truncatePolynomial(v, order_1, order_2), F):
    F_t := H_t:
    lprint("Current F", F):
    lprint("Current F_t", F_t):
  end if:
end do:

end proc:

testBasisConversion := proc (basis, order_1, order_2)
local output_1, output_2:
lprint("Basis conversion of ", basis, " from ", order_1, " to ", order_2):
output_1 := basisConversion(basis, order_1, order_2):
output_2 := Basis(basis, order_2):
lprint("Result obtained using implemented algorithm ", output_1):
lprint("Result obtained using Maple Basis algorithm ", output_2):
lprint("Are these results the same? ", evalb(output_1 = output_2)):
lprint(""):
end proc:

# ---------------------------------------------------------------------------
# Testing
#testBasisConversion([y^2-x,x^2-y*z-1,z^2-x], grlex(x, y, z), plex(x, y, z)):
#testBasisConversion([y^2-x,x^2-y*z-1,z^2-x], grlex(x, y, z), plex(x, z, y)):
#testBasisConversion([y^2-x,x^2-y*z-1,z^2-x], grlex(x, y, z), plex(y, x, z)):
#testBasisConversion([y^2-x,x^2-y*z-1,z^2-x], grlex(x, y, z), plex(y, z, x)):
#testBasisConversion([y^2-x,x^2-y*z-1,z^2-x], grlex(x, y, z), plex(z, x, y)):
#testBasisConversion([y^2-x,x^2-y*z-1,z^2-x], grlex(x, y, z), plex(z, y, x)):
#testBasisConversion([y^2-x,x^2-y*z-1,z^2-x], plex(x, y, z), grlex(x, y, z)):
#testBasisConversion([x*y+z-x*z, x^2-z, 2*x^3-x^2*y*z-1], tdeg(x, y, z), plex(z, y, x)):
#testBasisConversion([x + y + z, x*y + y*z + z*x, x*y*z-1], grlex(x, y, z), plex(x, y, z)):
testBasisConversion([x + y + z + u, x*y + y*z + z*u + u*x, x*y*z + y*z*u + z*u*x + u*x*y, x*y*z*u-1], grlex(x, y, z, u), plex(x, y, z, u)):
# ---------------------------------------------------------------------------
