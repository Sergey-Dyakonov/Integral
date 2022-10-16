@x_min = @y_min = @z_min = 0
@x_max = @y_max = @z_max = 1
@step = 0.0001
@a_param = 1
@n_param = 1
@k_param = 1

@func_1 = -> (a, x) { a * (1 - x) * x }
@func_2 = -> (n, y) { Math.exp(-n * y) }
@func_3 = -> (k, z) { Math.sin(Math::PI * k * z) }

direct_integration = ->{

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

def total_square(func)
  func.call(@x_min, @x_max, @step, @func_1, @a_param) *
    func.call(@y_min, @y_max, @step, @func_2, @n_param) *
    func.call(@z_min, @z_max, @step, @func_3, @k_param)
end

puts total_square(trapezium_method_integration)
puts total_square(square_method_integration)