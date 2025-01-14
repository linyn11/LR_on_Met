import tensorflow as tf

from absl.testing import parameterized
from tensorflow.python.keras import combinations

from tensorflow_riemopt.manifolds.test_invariants import TestInvariants
from tensorflow_riemopt.manifolds.euclidean import Euclidean


@combinations.generate(
    combinations.combine(
        mode=["graph", "eager"],
        manifold=[Euclidean()],
        shape=[(2,), (2, 2)],
        dtype=[tf.float32, tf.float64],
    )
)
class EuclideanTest(tf.test.TestCase, parameterized.TestCase):
    test_random = TestInvariants.check_random

    test_dist = TestInvariants.check_dist

    test_inner = TestInvariants.check_inner

    test_proj = TestInvariants.check_proj

    test_exp_log_inverse = TestInvariants.check_exp_log_inverse

    test_transp_retr = TestInvariants.check_transp_retr

    test_ptransp_inverse = TestInvariants.check_ptransp_inverse

    test_ptransp_inner = TestInvariants.check_ptransp_inner

    test_geodesic = TestInvariants.check_geodesic

    test_pairmean = TestInvariants.check_pairmean
