@x_min = @y_min = @z_min = 0
@x_max = @y_max = @z_max = 1
@step = 0.0001
@quantity = 100_000
@a_param = 1.0
@n_param = 1.0
@k_param = 1.0

@func_1 = -> x { @a_param.to_f * (1.0 - x.to_f) * x.to_f }
@func_2 = -> y { Math.exp(-@n_param.to_f * y.to_f) }
@func_3 = -> z { Math.sin(Math::PI * @k_param.to_f * z.to_f) }

analytical_eval = -> {
  func_1 = -> x { -@a_param * (x ** 2.0) * (2.0 * x - 3.0) / 6.0 }
  func_2 = -> y { -Math.exp(-@n_param * y) / @n_param }
  func_3 = -> z { -Math.cos(Math::PI * @k_param * z) / (Math::PI * @k_param) }
  (func_1.call(@x_min) - func_1.call(@x_max)).abs *
    (func_2.call(@y_min) - func_2.call(@y_max)).abs *
    (func_3.call(@z_min) - func_3.call(@z_max)).abs
}

trapezium_method_integration = ->(min, max, step, func) {
  sum = 0
  (min..max).step(step).each { |variable|
    sum += @step * 0.5 * (func.call(variable) + func.call(variable + @step))
  }
  sum
}

square_method_integration = ->(min, max, step, func) {
  sum = 0
  (min..max).step(step).each { |variable|
    sum += func.call(variable)
  }
  # number_of_nodes = (max - min) / step
  # (max - min) * sum / number_of_nodes = (max - min) * sum / (max - min) / step = sum / 1 / step = sum * step
  #
  # after basic arithmetic transformations the operations above turns into the operation below
  sum * step
}

the_simplest_monte_carlo = ->(min, max, quantity, func) {
  sum = 0
  (0...quantity).each {
    sum += func.call(rand(min.to_f...max.to_f))
  }
  (max - min) * sum / quantity
}

geometric_monte_carlo = ->(min, max, quantity, func) {
  min_max = find_min_and_max(min, max, quantity, func)
  f_min = min_max[0]
  f_max = min_max[1]
  n_1 = 0
  (0...quantity).each {
    x = min + (max - min) * rand(min.to_f..max.to_f)
    y = f_min + (f_max - f_min) * rand(min.to_f..max.to_f)
    if func.call(x) > y
      n_1 += 1
    end
  }
  print "Answer: " + ((max - min) * ((f_max - f_min) * n_1 / quantity + f_min)).to_s + "\n"
  (max - min) * ((f_max - f_min) * n_1 / quantity + f_min)
}

def find_min_and_max(min, max, quantity, func)
  f_min = f_max = func.call(min)
  (0...quantity).each { |val|
    f = func.call((max - min).to_f / quantity * val)
    if f < f_min
      f_min = f
    elsif f > f_max
      f_max = f
    end
  }
  [f_min, f_max]
end

def total_square(func)
  func.call(@x_min, @x_max, @step, @func_1) *
    func.call(@y_min, @y_max, @step, @func_2) *
    func.call(@z_min, @z_max, @step, @func_3)
end

def monte_carlo_total_square(func)
  func.call(@x_min, @x_max, @quantity, @func_1) *
    func.call(@y_min, @y_max, @quantity, @func_2) *
    func.call(@z_min, @z_max, @quantity, @func_3)
end

puts analytical_eval.call
puts total_square(trapezium_method_integration)
puts total_square(square_method_integration)
puts monte_carlo_total_square(the_simplest_monte_carlo)
puts monte_carlo_total_square(geometric_monte_carlo)