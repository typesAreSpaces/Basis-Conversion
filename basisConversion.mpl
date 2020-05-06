with(LinearAlgebra, 
  Multiply, Transpose):
with(Groebner,      
  Basis, InterReduce, LeadingTerm, LeadingMonomial, TestOrder): 

writeto("output.txt"):

truncatePolynomial := proc (poly, order_1, order_2)
local leading_coeff_1, leading_mon_1, leading_term, 
curr_index, polys:

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

basisConversion := proc (basis, order_1, order_2)
local F, F_t, num_iter, 
F_t_gb, multipliers, 
G, repeat, curr_index, num_elements:

lprint("Step 1"):
F   := Basis(basis, order_1):
F_t := map(v -> truncatePolynomial(v, order_1, order_2), F):
lprint("Current F", F):
lprint("Current F_t", F_t):

num_iter := 1:

while true do
  lprint("Iteration ", num_iter):
  num_iter := num_iter + 1:

  lprint("Step 3 & 4"):
  F_t_gb, multipliers := Basis(F_t, order_2, output=extended):
  lprint("Matrix of multipliers M': ", multipliers):
  lprint("H (= M' F_t): ", F_t_gb):
  
  lprint("Step 5"):
  G := convert(simplify(
  Multiply(convert(multipliers, Matrix), Transpose(convert(F, Matrix)))), list):
  lprint("G (= M' * F): ", G):

  lprint("Step 6 and 7"):
  repeat       := false:
  curr_index   := 1:
  num_elements := numelems(F_t_gb):
  while curr_index < num_elements do
    if evalb(LeadingMonomial(F_t_gb[curr_index], order_2) 
      <> LeadingMonomial(G[curr_index], order_2)) then
      lprint("There are the witness polynomials that prove"\
       "G_s is not empty (g_j, h_j) respectively ", 
      F_t_gb[curr_index], G[curr_index]):
      repeat := true:
      break:
    end if:
    curr_index := curr_index + 1:
  end do:

  #F   := G:
  F   := InterReduce(G, order_2):
  if repeat = false then
    lprint("Done"):
    return F:
  else
    F_t := map(v -> truncatePolynomial(v, order_1, order_2), F):
    lprint("Current F", F):
    lprint("Current F_t", F_t):
  end if:
end do:

end proc:

testBasisConversion := proc (basis, order_1, order_2)
lprint("Basis conversion of ", basis, " from ", order_1, " to ", order_2):
lprint("Result obtained using implemented algorithm ", basisConversion(basis, order_1, order_2)):
lprint("Result obtained using Maple Basis algorithm ", Basis(basis, order_2)):
lprint(""):
end proc:


# ---------------------------------------------------------------------------
# Test
testBasisConversion([y^2-x,x^2-y*z-1,z^2-x], grlex(x, y, z), plex(x, y, z)):
testBasisConversion([y^2-x,x^2-y*z-1,z^2-x], plex(x, y, z), grlex(x, y, z)):
# ---------------------------------------------------------------------------
