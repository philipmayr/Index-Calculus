# p is a safe prime (two times a prime number plus one)
# p = 1283
p = 232487
# g is a primitive root modulo p
# g = 264
g = 229401
# h is a prime modulated power of primitive root g
# h = 910
h = 171020

# exponential index
x = nil

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

smooth_number_logarithms.each do |relation|
    puts "i æ #{relation[:index]}, number æ #{relation[:number]}, exponents æ #{relation[:exponent_multiplicity_counts].inspect}"
end
