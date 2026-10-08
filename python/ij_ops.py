"""Python equivalents of the ImageJ operations the loop-calling macro uses.

Each function follows the ImageJ 1.54 source (ij/plugin/filter/*.java, ij/process/*.java) so
that bactoloops.py reproduces the macro as closely as possible. Images are float32 numpy arrays
indexed [row, column] (ImageJ's [y, x]).
"""
import numpy as np
from scipy import ndimage as ndi


# ---------------------------------------------------------------- Subtract Background
def _filter3x3_mean(a):
    """BackgroundSubtracter.filter3x3(MEAN): 3-pixel mean along rows, then along columns,
    edge pixels repeated, in float32 as ImageJ does."""
    a = a.astype(np.float32)
    for axis in (1, 0):
        p = np.pad(a, [(1, 1) if ax == axis else (0, 0) for ax in (0, 1)], mode="edge")
        if axis == 1:
            a = (p[:, :-2] + p[:, 1:-1] + p[:, 2:]) * np.float32(0.333333333)
        else:
            a = (p[:-2] + p[1:-1] + p[2:]) * np.float32(0.333333333)
    return a.astype(np.float32)


def _ball(radius):
    """RollingBall: shrink factor and the (trimmed) ball patch."""
    if radius <= 10:
        shrink, trim = 1, 24
    elif radius <= 30:
        shrink, trim = 2, 24
    elif radius <= 100:
        shrink, trim = 4, 32
    else:
        shrink, trim = 8, 40
    r = max(radius / shrink, 1)
    xtrim = int(trim * r) // 100
    half = int(round(r - xtrim))
    y, x = np.mgrid[-half:half + 1, -half:half + 1]
    t = r * r - x * x - y * y
    return shrink, np.where(t > 0, np.sqrt(np.maximum(t, 0)), 0).astype(np.float32)


def _interp_arrays(length, small_length, shrink):
    """BackgroundSubtracter.makeInterpolationArrays (Java integer division truncates toward 0)."""
    i = np.arange(length)
    idx = np.trunc((i - shrink // 2) / shrink).astype(int)
    idx = np.minimum(idx, small_length - 2)
    dist = (i + 0.5) / shrink - (idx + 0.5)
    return idx, (1 - dist).astype(np.float32)


def rolling_ball_background(img, radius, presmooth=True):
    """BackgroundSubtracter.rollingBallFloatBackground (no inversion): the background image."""
    a = img.astype(np.float32)
    if presmooth:
        a = _filter3x3_mean(a)
    shrink, ball = _ball(radius)
    h, w = a.shape
    if shrink > 1:
        sh, sw = -(-h // shrink), -(-w // shrink)
        p = np.full((sh * shrink, sw * shrink), np.inf, np.float32)
        p[:h, :w] = a
        small = p.reshape(sh, shrink, sw, shrink).min(axis=(1, 3))
    else:
        small = a
    # roll the ball: grey opening with the ball patch; ball centres may sit up to `half` outside
    # the image, and only in-image points count when deciding how high the ball can rise
    half = ball.shape[0] // 2
    padded = np.pad(small, half, mode="constant", constant_values=np.inf)
    z = ndi.grey_erosion(padded, structure=ball, mode="constant", cval=np.inf)
    z[~np.isfinite(z)] = -np.inf
    bg = ndi.grey_dilation(z, structure=ball, mode="constant", cval=-np.inf)[half:-half, half:-half]
    bg = bg.astype(np.float32)
    if shrink > 1:
        sh, sw = bg.shape
        yi, yw = _interp_arrays(h, sh, shrink)
        xi, xw = _interp_arrays(w, sw, shrink)
        rows = bg[yi][:, xi] * xw[None, :] + bg[yi][:, xi + 1] * (1 - xw[None, :])
        rows1 = bg[yi + 1][:, xi] * xw[None, :] + bg[yi + 1][:, xi + 1] * (1 - xw[None, :])
        bg = (rows * yw[:, None] + rows1 * (1 - yw[:, None])).astype(np.float32)
    return bg


def _translate(a, dx, dy):
    """Image > Transform > Translate, interpolation None: shift by (dx, dy) pixels, vacated pixels 0."""
    out = np.zeros_like(a)
    h, w = a.shape
    ys, yd = (slice(0, h - dy), slice(dy, h)) if dy >= 0 else (slice(-dy, h), slice(0, h + dy))
    xs, xd = (slice(0, w - dx), slice(dx, w)) if dx >= 0 else (slice(-dx, w), slice(0, w + dx))
    out[yd, xd] = a[ys, xs]
    return out


def periodic_background(padded, n, pad_half, radius, phases=4, margin=100):
    """The macro's rotation-invariant local background (step 7c, localBackgroundPhaseAverage).

    ImageJ's rolling ball shrinks the image on a grid anchored at pixel (0, 0), so its background
    depends on where a circular map starts. The map core (n x n, at pad_half in `padded`) is tiled
    3 x 3, cropped to the padded frame plus `margin` px, and the rolling-ball background is averaged
    over every (x, y) grid phase, then made exactly symmetric (mean with its transpose). Returns the
    background in the padded frame, in float32 as ImageJ computes it.
    """
    core = padded[pad_half:pad_half + n, pad_half:pad_half + n].astype(np.float32)
    size = padded.shape[0] + 2 * margin
    o = n - pad_half - margin
    tile = np.tile(core, (3, 3))[o:o + size, o:o + size]
    acc = np.zeros((size, size), np.float32)
    for phy in range(phases):
        for phx in range(phases):
            t = _translate(tile, phx, phy)
            b = rolling_ball_background(t, radius)
            acc = (acc + _translate(b, -phx, -phy)).astype(np.float32)
    acc = (acc / np.float32(phases * phases)).astype(np.float32)
    bg = acc[margin:margin + padded.shape[0], margin:margin + padded.shape[1]]
    return ((bg + bg.T) / np.float32(2)).astype(np.float32)


def circular_subtract_background(core, radius, pad=100):
    """Rotation-invariant Subtract Background for an n x n circular map (no padding needed): the map
    minus periodic_background (all grid phases of ImageJ's shrink factor for this radius, symmetrized).
    Use this instead of subtract_background(np.pad(core, pad, mode="wrap"), radius)[pad:-pad, pad:-pad]."""
    core = np.asarray(core, np.float32)
    n = core.shape[0]
    padded = np.pad(core, pad, mode="wrap")
    bg = periodic_background(padded, n, pad, radius, phases=_ball(radius)[0])
    return (padded - bg)[pad:pad + n, pad:pad + n].astype(np.float32)


def subtract_background(img, radius):
    """Process > Subtract Background, rolling ball, default options, on a 32-bit image."""
    return (img.astype(np.float32) - rolling_ball_background(img, radius)).astype(np.float32)


# ---------------------------------------------------------------- rank and Gaussian filters
def median(img, radius):
    """Process > Filters > Median (RankFilters): circular kernel of ImageJ's shape, edge pixels
    repeated outside the image."""
    r = radius
    if 1.5 <= r < 1.75:
        r = 1.75
    elif 2.5 <= r < 2.85:
        r = 2.85
    r2 = int(r * r) + 1
    k = int(np.sqrt(r2 + 1e-10))
    fp = np.zeros((2 * k + 1, 2 * k + 1), bool)
    fp[k, :] = True
    for y in range(1, k + 1):
        dx = int(np.sqrt(r2 - y * y + 1e-10))
        fp[k - y, k - dx:k + dx + 1] = True
        fp[k + y, k - dx:k + dx + 1] = True
    return ndi.median_filter(img.astype(np.float32), footprint=fp, mode="nearest")


def _gauss_kernel(sigma, accuracy, max_radius=50):
    """GaussianBlur.makeGaussianKernel: one-sided kernel with ImageJ's edge correction."""
    kr = int(np.ceil(sigma * np.sqrt(-2 * np.log(accuracy)))) + 1
    max_radius = max(max_radius, 50)
    kr = min(kr, max_radius)
    kern = np.exp(-0.5 * np.arange(kr) ** 2 / sigma / sigma).astype(np.float32).astype(float)
    if kr < max_radius and kr > 3:
        sqrt_slope = np.inf
        r = kr
        while r > kr // 2:
            r -= 1
            a = np.sqrt(kern[r]) / (kr - r)
            if a < sqrt_slope:
                sqrt_slope = a
            else:
                break
        for r1 in range(r + 2, kr):
            kern[r1] = np.float32((kr - r1) ** 2 * sqrt_slope * sqrt_slope)
    total = kern[0] + 2 * kern[1:].sum() if kr < max_radius else sigma * np.sqrt(2 * np.pi)
    return (kern / total).astype(np.float32)


def gaussian(img, sigma, bit_depth=32):
    """Process > Filters > Gaussian Blur (sigma in pixels): x then y, edge pixels repeated.
    8-bit/16-bit results are rounded back to integers as ImageJ does."""
    acc = 0.002 if bit_depth == 8 else 0.0002
    k = _gauss_kernel(sigma, acc, max(img.shape))
    full = np.r_[k[:0:-1], k].astype(np.float32)
    out = img.astype(np.float32)
    out = ndi.correlate1d(out, full, axis=1, mode="nearest")
    out = ndi.correlate1d(out, full, axis=0, mode="nearest")
    if bit_depth == 16:
        out = np.clip(np.floor(out + 0.5), 0, 65535)
    elif bit_depth == 8:
        out = np.clip(np.floor(out + 0.5), 0, 255)
    return out.astype(np.float32)


# ---------------------------------------------------------------- Image > Adjust > Size (bilinear)
def resize_bilinear_16bit(img, factor):
    """ShortProcessor.resize with BILINEAR interpolation, as run("Size...", ...) does: the
    destination centre is shifted by scale/4 (the source of ImageJ's ~0.14-pixel offset), and
    values are rounded to 16-bit."""
    h, w = img.shape
    H, W = h * factor, w * factor
    def coords(n_dst, n_src):
        scale = n_dst / n_src
        c = n_dst / 2.0 + scale / 4.0
        s = (np.arange(n_dst) - c) / scale + n_src / 2.0
        s = np.where(s < 0, 0.0, s)
        s = np.where(s >= n_src - 1.0, n_src - 1.001, s)
        return s
    ys, xs = coords(H, h), coords(W, w)
    yb, xb = ys.astype(int), xs.astype(int)
    yf, xf = (ys - yb)[:, None], (xs - xb)[None, :]
    a = img.astype(float)
    ll = a[np.ix_(yb, xb)]; lr = a[np.ix_(yb, xb + 1)]
    ul = a[np.ix_(yb + 1, xb)]; ur = a[np.ix_(yb + 1, xb + 1)]
    upper = ul + xf * (ur - ul)
    lower = ll + xf * (lr - ll)
    v = lower + yf * (upper - lower)
    return np.clip(np.floor(v + 0.5), 0, 65535).astype(np.float32)


# ---------------------------------------------------------------- binary and particle tools
def histogram_quantiles(img, qs, n_bins=65536):
    """The macro's computeRobustStats/getHistogram: quantiles read off an n_bins histogram
    spanning the image's own [min, max] (bin lower edges, as the macro reports them)."""
    v = img.ravel().astype(float)
    lo, hi = v.min(), v.max()
    if hi <= lo:
        return [lo for _ in qs]
    bw = (hi - lo) / n_bins
    idx = np.minimum(((v - lo) / bw).astype(int), n_bins - 1)
    cum = np.cumsum(np.bincount(idx, minlength=n_bins))
    return [lo + np.searchsorted(cum, q * len(v)) * bw for q in qs]


def particles(binary, min_area=0, max_area=np.inf, area_scale=1.0):
    """Analyze Particles on a binary image: 8-connected particles, in the raster order ImageJ
    finds them (first pixel top-to-bottom, left-to-right), with pixel-count areas (x area_scale)
    between min_area and max_area. Returns a label image and the kept labels in that order."""
    lab, n = ndi.label(binary, structure=np.ones((3, 3), bool))
    if n == 0:
        return lab, []
    first = ndi.minimum(np.arange(lab.size).reshape(lab.shape), lab, index=np.arange(1, n + 1))
    area = ndi.sum(np.ones_like(lab), lab, index=np.arange(1, n + 1)) * area_scale
    order = np.argsort(first)
    keep = [int(i + 1) for i in order if min_area <= area[i] <= max_area]
    return lab, keep


def watershed_split(binary):
    """Process > Binary > Watershed: ImageJ's EDM + MaximumFinder segmentation, ported exactly
    in ij_watershed.py (needs numba)."""
    from ij_watershed import watershed
    return watershed(binary)
