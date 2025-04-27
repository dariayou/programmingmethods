import java.util.*;
import java.util.stream.Collectors; 

class Root {
    double value;

    public Root(double value) { 
        this.value = value;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Root root = (Root) o;
        return Double.compare(root.value, value) == 0;
    }

    @Override
    public int hashCode() { 
        long temp = Double.doubleToLongBits(value);
        return (int) (temp ^ (temp >>> 32)); 
    }
}

public class Main {

    static double calc(double x, double[] c) { 
        double r = 0;
        for (int i = 0; i < c.length; i++) {
            r += c[i] * Math.pow(x, c.length - 1 - i);
        }
        return r;
    }

    static double find(double a, double b, double[] c, double p) { 
        while (Math.abs(b - a) > p) {
            double m = (a + b) / 2.0;
            if (calc(a, c) * calc(m, c) <= 0) {
                b = m;
            } else {
                a = m;
            }
        }
        return (a + b) / 2.0;
    }

    static List<Root> solve(double[] r, double p, double[] c) {  
        List<Root> roots = new ArrayList<>();
        double s = (r[1] - r[0]) / 100.0;

        for (double x = r[0]; x <= r[1]; x += s) {
            double n = x + s;
            if (calc(x, c) * calc(n, c) <= 0) {
                double root = find(x, n, c, p);
                int d = (int) -Math.log10(p);
                root = Math.round(root * Math.pow(10, d)) / Math.pow(10, d);

                roots.add(new Root(root));
            }
        }

        Set<Root> unique = new HashSet<>(roots); 
        roots.clear();
        roots.addAll(unique);

        return roots;
    }

    public static void main(String[] args) {
        double[] r1 = {-5, 5};
        double p1 = 0.0001;
        double[] c1 = {1, -1, -6};
        List<Root> res1 = solve(r1, p1, c1);
        System.out.println("Пример 1: y = x^2 - x - 6");
        System.out.println("Найденные корни: " + res1.stream().map(root -> String.valueOf(root.value)).collect(Collectors.toList()));
        System.out.println("Ожидаемый результат: [-2.0, 3.0]");
        System.out.println();

        double[] r2 = {-5, 5};
        double p2 = 0.0001;
        double[] c2 = {1, 0, 0, 0};
        List<Root> res2 = solve(r2, p2, c2);
        System.out.println("Пример 2: y = x^3");
        System.out.println("Найденные корни: " + res2.stream().map(root -> String.valueOf(root.value)).collect(Collectors.toList()));
        System.out.println("Ожидаемый результат: [0.0]");
    }
}
