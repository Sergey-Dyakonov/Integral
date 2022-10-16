@x_min = @y_min = @z_min = 0
@x_max = @y_max = @z_max = 1
@step = 0.0001
@quantity = 100_000
@a_param = 1
@n_param = 1
@k_param = 1

@func_1 = -> (a, x) { a * (1 - x) * x }
@func_2 = -> (n, y) { Math.exp(-n * y) }
@func_3 = -> (k, z) { Math.sin(Math::PI * k * z) }

analytical_eval = -> {
  func_1 = -> (a, x) { -a * (x ** 2.0) * (2.0 * x - 3.0) / 6.0 }
  func_2 = -> (b, y) { -Math.exp(-b * y) / b }
  func_3 = -> (c, z) { -Math.cos(Math::PI * c * z) / (Math::PI * c) }
  (func_1.call(@a_param, @x_min) - func_1.call(@a_param, @x_max)).abs *
    (func_2.call(@n_param, @y_min) - func_2.call(@n_param, @y_max)).abs *
    (func_3.call(@k_param, @z_min) - func_3.call(@k_param, @z_max)).abs
}

trapezium_method_integration = ->(min, max, step, func, param) {
  sum = 0
  (min + step...max).step(step).each { |variable|
    sum += @step * 0.5 * (func.call(param, variable) + func.call(param, variable + @step))
  }
  sum
}

square_method_integration = ->(min, max, step, func, param) {
  sum = 0
  (min + step...max).step(step).each { |variable|
    sum += func.call(param, variable)
  }
  # number_of_nodes = (max - min) / step
  # (max - min) * sum / number_of_nodes = (max - min) * sum / (max - min) / step = sum / 1 / step = sum * step
  #
  # after basic arithmetic transformations the operations above turns into the operation below
  sum * step
}

the_simplest_monte_carlo = ->(min, max, quantity, func, param) {
  sum =0
  (0...quantity).each {
    sum += func.call(param, rand(min.to_f...max.to_f))
  }
  (max - min) * sum / quantity
}

def total_square(func)
  func.call(@x_min, @x_max, @step, @func_1, @a_param) *
    func.call(@y_min, @y_max, @step, @func_2, @n_param) *
    func.call(@z_min, @z_max, @step, @func_3, @k_param)
end

def monte_carlo_total_square(func)
  func.call(@x_min, @x_max, @quantity, @func_1, @a_param) *
    func.call(@y_min, @y_max, @quantity, @func_2, @n_param) *
    func.call(@z_min, @z_max, @quantity, @func_3, @k_param)
end

puts analytical_eval.call
puts total_square(trapezium_method_integration)
puts total_square(square_method_integration)
puts monte_carlo_total_square(the_simplest_monte_carlo)