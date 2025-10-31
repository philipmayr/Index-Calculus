# p is a safe prime (two times a prime number plus one)
# p = 1283
p = 232487
# g is a primitive root modulo p: a primitive root is a generator of all elements in the group
# all elements can be written as powers of a primitive root (generator) modulo p
# a number g is a primitive root modulo p if g is coprime to p
# the powers of g modulo p generate all elements of the multiplicative group modulo p
# the multiplicative group modulo p is the set of integers from 1 to p - 1 that are invertible under multiplication modulo p
# an element is invertible modulo p if it is coprime to p
# g = 264
g = 229401
# h is a prime modulated power of primitive root g
# h = 910
h = 171020

# exponential index
x = 0

# factor base 𝓕
factor_base = [2, 3, 5, 7]

smooth_number_logarithms = Array.new()

def exponentiate_modularly(base, exponent, modulus)
    return 0 if base == 0 || modulus == 1 
    return 1 if exponent == 0
        
    base %= modulus
    
    power = 1
    
    while exponent > 0
        power = (power * base) % modulus if exponent & 1 == 1
        
        base = (base * base) % modulus
        exponent >>= 1
    end
    
    return power
end

def factor_number_over_base(number, factor_base)
    leftover_cofactor = number
    exponent_multiplicity_counts = Array.new(factor_base.length, 0)
    
    return [false, exponent_multiplicity_counts, leftover_cofactor] if leftover_cofactor == 0
    return [true, exponent_multiplicity_counts, 1] if leftover_cofactor == 1
    
    factor_base.each_with_index do |prime, index|
        while leftover_cofactor % prime == 0
            leftover_cofactor /= prime
            exponent_multiplicity_counts[index] += 1
        end
    end
    
    [leftover_cofactor == 1, exponent_multiplicity_counts, leftover_cofactor]
end

for i in 1..1000
    modular_power = exponentiate_modularly(g, i, p)
    is_smooth, exponent_multiplicity_counts, leftover_cofactor = factor_number_over_base(modular_power, factor_base)
    
    if is_smooth
        smooth_number_logarithms << { index: i, number: modular_power, exponent_multiplicity_counts: exponent_multiplicity_counts }
    end
end

EMC = exponent_multiplicity_counts

# puts EMC

# smooth_number_logarithms.each do |relation|
#     puts "i æ #{relation[:index]}, number æ #{relation[:number]}, exponents æ #{relation[:exponent_multiplicity_counts].inspect}"
# end

# find j where gʲh is factorable over 𝓕

for j in 1..p - 1
    modular_power_of_g = exponentiate_modularly(g, j, p)
    modular_power_of_g_multiple_of_h = (modular_power_of_g * h) % p
    
    is_smooth, exponent_multiplicity_counts, leftover_cofactor = factor_number_over_base(modular_power_of_g_multiple_of_h, factor_base)
    
    if is_smooth
        break
    end
end

modulus = 116243

matrix = [[5, 2,  0,  3,  482],
          [2, 9,  0,  0,  600],
          [9, 4,  0,  0,  876],
          [6, 0,  2,  0,  948]]
          
matrix[0] = smooth_number_logarithms[0][:exponent_multiplicity_counts].dup.append(smooth_number_logarithms[0][:index])
matrix[1] = smooth_number_logarithms[1][:exponent_multiplicity_counts].dup.append(smooth_number_logarithms[1][:index])
matrix[2] = smooth_number_logarithms[2][:exponent_multiplicity_counts].dup.append(smooth_number_logarithms[2][:index])
matrix[3] = smooth_number_logarithms[3][:exponent_multiplicity_counts].dup.append(smooth_number_logarithms[3][:index])

# print matrix[0]

def modulate(number, modulus)
    number %= modulus
    number += modulus if number < 0
    number
end

def find_modular_multiplicative_inverse(base, modulus)
    base %= modulus
    raise "No inverse for 0 modulo #{modulus}" if base == 0
    t0, t1 = 0, 1
    r0, r1 = modulus, base
    while r1 != 0
        q = r0 / r1
        r0, r1 = r1, r0 - q * r1
        t0, t1 = t1, t0 - q * t1
    end
    t0 % modulus
end

# ERO1
def swap_roles(first_row, second_row)
    return second_row.dup, first_row.dup
end

# ERO2
def scale_row(row, constant)
    row.map { |element| modulate(element.to_i * constant.to_i, modulus) }
end

# ERO3
def add_multiplied_row(source_row, target_row, constant, modulus)
    target_row.each_with_index.map do |element, index|
        modulate(element.to_i + constant.to_i * source_row[index].to_i, modulus)
    end
end

def ensure_pivot_found_in_matrix!(matrix, column, modulus)
    if matrix[column][column] % modulus == 0
        swap_row_index = (column + 1...matrix.length).find { |row_index| matrix[row_index][column] % modulus != 0 }
        if swap_row_index
            matrix[column], matrix[swap_row_index] = swap_roles(matrix[column], matrix[swap_row_index])
        else
            raise "Zero pivot found in column #{column} modulo #{modulus}."
        end
    end
end

matrix.map! { |row| row.map { |element| modulate(element, modulus) } }

# --- Forward elimination
n = matrix.length
(0...n - 1).each do |column|
    ensure_pivot_found_in_matrix!(matrix, column, modulus)
    pivot = matrix[column][column] % modulus
    inverse_pivot = find_modular_multiplicative_inverse(pivot, modulus)
    ((column + 1)...n).each do |row_index|
        # multiplier = matrix[row_index][column] / pivot  (mod modulus)
        multiplier = (matrix[row_index][column] * inverse_pivot) % modulus
        # new_row = target_row - multiplier * pivot_row  (mod modulus)
        matrix[row_index] = add_multiplied_row(matrix[column], matrix[row_index], (-multiplier) % modulus, modulus)
    end
end

# --- Backward substitution
solution = Array.new(n, 0)
(n - 1).downto(0) do |i|
    right_hand_side = matrix[i][n] % modulus
    sum_of_products = 0
    ((i + 1)...n).each { |j| sum_of_products = (sum_of_products + matrix[i][j] * solution[j]) % modulus }
    right_hand_side_minus_sum = (right_hand_side - sum_of_products) % modulus
    inverse_diagonal = find_modular_multiplicative_inverse(matrix[i][i] % modulus, modulus)
    solution[i] = (right_hand_side_minus_sum * inverse_diagonal) % modulus
end

# print matrix and solutions
# print matrix[0]
# puts
# print matrix[1]
# puts
# print matrix[2]
# puts
# print matrix[3]
# puts

# print solution

congruences = solution.dup

for i in 0...exponent_multiplicity_counts.size
    x += congruences[i] * exponent_multiplicity_counts[i]
end

x -= j

puts x % (p - 1)
