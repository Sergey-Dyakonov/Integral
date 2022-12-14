require 'benchmark'

@x_min = @y_min = @z_min = 0
@x_max = @y_max = @z_max = 1
@step = 0.0001
@quantity = 100_000
@a_param = @n_param = @k_param = 1.0

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
  (min..max).step(step) { |variable|
    sum += step * 0.5 * (func.call(variable) + func.call(variable + step))
  }
  sum
}

square_method_integration = ->(min, max, step, func) {
  sum = 0
  (min..max).step(step) { |variable|
    sum += func.call(variable)
  }
  # (max - min) * sum / number_of_nodes = (max - min) * sum / (max - min) / step = sum / 1 / step = sum * step
  sum * step
}

runge = -> {
  mistake = 0
  s_h = square_method_integration.call(@x_min, @x_max, @step, @func_1)
  s_2h = square_method_integration.call(@x_min, @x_max, 2 * @step, @func_1)
  mistake += s_h - s_2h
  s_h = square_method_integration.call(@y_min, @y_max, @step, @func_2)
  s_2h = square_method_integration.call(@y_min, @y_max, 2 * @step, @func_2)
  mistake += s_h - s_2h
  s_h = square_method_integration.call(@z_min, @z_max, @step, @func_3)
  s_2h = square_method_integration.call(@z_min, @z_max, 2 * @step, @func_3)
  mistake += s_h - s_2h
  (mistake / 3.0).abs
}

the_simplest_monte_carlo = ->(min, max, func) {
  sum = 0
  (0...@quantity).each {
    sum += func.call(rand(min.to_f...max.to_f))
  }
  (max - min) * sum / @quantity
}

@simple_dispersion = ->(min, max, func) {
  sum = 0
  sum_2 = 0
  (0...@quantity).each {
    f = func.call(rand(min.to_f...max.to_f))
    sum += f
    sum_2 += f ** 2
  }
  (sum_2 / @quantity - (sum / @quantity) ** 2) / @quantity
}

def simple_sigma
  sigma_1 = (@x_max - @x_min) * Math.sqrt(@simple_dispersion.call(@x_min, @x_max, @func_1))
  sigma_2 = (@y_max - @y_min) * Math.sqrt(@simple_dispersion.call(@y_min, @y_max, @func_2))
  sigma_3 = (@z_max - @z_min) * Math.sqrt(@simple_dispersion.call(@z_min, @z_max, @func_3))
  (sigma_1 + sigma_2 + sigma_3) / 3.0
end

geometric_monte_carlo = ->(min, max, func) {
  min_max = find_min_and_max(min, max, func)
  f_min = min_max[0]
  f_max = min_max[1]
  n_1 = 0
  (0...@quantity).each {
    x = min + (max - min) * rand(min.to_f..max.to_f)
    y = f_min + (f_max - f_min) * rand(min.to_f..max.to_f)
    if func.call(x) > y
      n_1 += 1
    end
  }
  (max - min) * ((f_max - f_min) * n_1 / @quantity + f_min)
}

epsilon = -> (min, max, func) {
  x_i = rand(min.to_f..max.to_f)
  y_i = rand(min.to_f..max.to_f)
  func.call(x_i) > y_i ? 0 : 1
}

@geom_dispersion = ->(min, max, func) {
  sum = 0
  (0...@quantity).each {
    sum += epsilon.call(min, max, func)
  }
  sum.to_f / @quantity * (1.0 - sum.to_f / @quantity)
}

def geometric_sigma
  min_max = find_min_and_max(@x_min, @x_max, @func_1)
  sigma_1 = (@x_max - @x_min) * (min_max[1] - min_max[0]) * Math.sqrt(@geom_dispersion.call(@x_min, @x_max, @func_1) / @quantity)
  min_max = find_min_and_max(@y_min, @y_max, @func_2)
  sigma_2 = (@y_max - @y_min) * (min_max[1] - min_max[0]) * Math.sqrt(@geom_dispersion.call(@y_min, @y_max, @func_2) / @quantity)
  min_max = find_min_and_max(@z_min, @z_max, @func_3)
  sigma_3 = (@z_max - @z_min) * (min_max[1] - min_max[0]) * Math.sqrt(@geom_dispersion.call(@z_min, @z_max, @func_3) / @quantity)
  (sigma_1 + sigma_2 + sigma_3) / 3.0
end

def find_min_and_max(min, max, func)
  f_min = f_max = func.call(min)
  (0...@quantity).each { |val|
    f = func.call((max - min).to_f / @quantity * val)
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
  func.call(@x_min, @x_max, @func_1) *
    func.call(@y_min, @y_max, @func_2) *
    func.call(@z_min, @z_max, @func_3)
end

puts "?????????? ???????????????????????? ????????????????: \n
?????????????? ????????????????????????: #{@x_min} < ?? < #{@x_max}; #{@y_min} < y < #{@y_max}; #{@z_min} < z < #{@z_max}
???????? ?????? ???????????????????????????????? ??????????????: #{@step}
?????????????????? ?????????????????? ?????? ???????????????????????? ??????????????: #{@quantity}
???????????????????????????? ??????????????:\tf1(x) = a * (1-x)x
\t\t\tf2(y) = exp(-ny)
\t\t\tf3(z) = sin(PIkz)
\t\t\tF(x,y,z) = f1(x) * f2(y) * f3(z)\n
?????????????? ??????????????????"
print "a ="
@a_param = gets.chomp.to_f
print "n ="
@n_param = gets.chomp.to_f
print "k ="
@k_param = gets.chomp.to_f

analytical_res = analytical_eval.call
square_res = total_square(square_method_integration)
square_bench = Benchmark.measure { total_square(square_method_integration) }
simple_monte_res = monte_carlo_total_square(the_simplest_monte_carlo)
sim_monte_bench = Benchmark.measure { monte_carlo_total_square(the_simplest_monte_carlo) }
geometric_monte_res = monte_carlo_total_square(geometric_monte_carlo)
geom_monte_bench = Benchmark.measure { monte_carlo_total_square(geometric_monte_carlo) }

puts "?????????????????????? ??????????: \t\t\t#{analytical_res}"
puts "?????????? ????????????????: \t\t\t#{total_square(trapezium_method_integration)}"
puts "?????????? ?????????????????????????? (1): \t\t#{square_res}\t??????????????: #{(analytical_res - square_res).abs}\t??????????????: #{runge.call}"
puts "???????????????????????? ?????????? ?????????? ?????????? (2): \t#{simple_monte_res}\t??????????????: #{(analytical_res - simple_monte_res).abs}\t??????????????: #{simple_sigma}"
puts "???????????????????????? ?????????? ?????????? ?????????? (3): \t#{geometric_monte_res}\t??????????????: #{(analytical_res - geometric_monte_res).abs}\t??????????????: #{geometric_sigma}"
puts "?????? ???????????????????? (1): #{square_bench.real} c"
puts "?????? ???????????????????? (2): #{sim_monte_bench.real} c"
puts "?????? ???????????????????? (3): #{geom_monte_bench.real} c"