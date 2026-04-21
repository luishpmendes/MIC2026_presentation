import math
import matplotlib.pyplot as plt
import seaborn as sns
import random
import sys
import statistics


def polynomial_mutation(allele: float, mutation_distribution: float) -> float:
    y = allele
    inner_exponent = mutation_distribution + 1.0
    outer_exponent = 1.0 / (mutation_distribution + 1.0)
    delta_l = y - 0.0
    delta_r = 1.0 - y
    delta = 0.0
    u = random.random()

    if u < 0.5:
        delta = math.pow(2.0 * u + (1.0 - 2.0 * u) * math.pow(1 - delta_l, inner_exponent), outer_exponent) - 1.0
    else:
        delta = 1.0 - math.pow(2.0 * (1.0 - u) + 2.0 * (u - 0.5) * math.pow(1 - delta_r, inner_exponent), outer_exponent)

    allele_prime = allele + delta

    if allele_prime < 0.0:
        allele_prime = 0.0
    elif allele_prime >= 1.0:
        allele_prime = 1.0 - 2 * sys.float_info.epsilon

    return allele_prime


def show(x: float = 0.38, mutation_distribution: float = 50) -> None:
    X = [x for _ in range(10000)]
    X_prime = [polynomial_mutation(x, mutation_distribution) for x in X]
    #print mean and std deviation
    print("X_prime")
    print("Mean: ", statistics.mean(X_prime))
    print("Std Dev: ", statistics.stdev(X_prime))
    print("Max: ", max(X_prime))
    print("Min: ", min(X_prime))
    print("Median: ", statistics.median(X_prime))
    print("Mode: ", statistics.mode(X_prime))
    print("1st Quartile: ", statistics.quantiles(X_prime, n=4)[0])
    print("3rd Quartile: ", statistics.quantiles(X_prime, n=4)[2])
    X_delta = [abs(X_prime[i] - X[i]) for i in range(len(X))]
    print("X_delta")
    print("Mean: ", statistics.mean(X_delta))
    print("Std Dev: ", statistics.stdev(X_delta))
    print("Max: ", max(X_delta))
    print("Min: ", min(X_delta))
    print("Median: ", statistics.median(X_delta))
    print("Mode: ", statistics.mode(X_delta))
    print("1st Quartile: ", statistics.quantiles(X_delta, n=4)[0])
    print("3rd Quartile: ", statistics.quantiles(X_delta, n=4)[2])
    plt.figure(figsize=(7, 5))
    plt.title(fontsize="x-large", label="Polynomial Mutation")
    plt.xlabel(fontsize="x-large", xlabel="Allele Value")
    plt.ylabel(fontsize="x-large", ylabel="Probability")
    plt.xlim(0, 1)  # Set x range from 0 to 1
    plt.tick_params(axis="both", which="both", labelsize="x-large")
    sns.histplot(X_prime, color="#7F72B4", label="Mutated Distribution", stat="probability", linewidth=0, kde=True)
    plt.grid(alpha=0.5, color='gray', linestyle='dashed', linewidth=0.5, which='both')
    plt.axvline(x=x, color='#E5761A', linestyle='--', label="Original Value", linewidth=2)
    plt.legend(fontsize="x-large")
    plt.tight_layout()
    plt.savefig("figs/polynomial_mutation.pdf", format="pdf")
    plt.savefig("figs/polynomial_mutation.png", format="png")
    plt.show()


if __name__ == "__main__":
    seed = 3429721499
    random.seed(seed)
    discard_count = 60 * 60
    show()
